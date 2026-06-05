import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report

# 1. قراءة البيانات من ملفات epitopes
files = {
    "ELAGIGILTV": "VDJdb_ELAGIGILTV_TRA_TRB_Input.csv",
    "AAGIGILTV": "VDJdb_AAGIGILTV_TRA_TRB_Input.csv",
    "IMDQVPFSV": "VDJdb_IMDQVPFSV_TRA_TRB_Input.csv",
    "KTWGQYWQV": "VDJdb_KTWGQYWQV_TRA_TRB_Input.csv",
    "YLEPGPVTA": "VDJdb_YLEPGPVTA_TRA_TRB_Input.csv"
}

dfs = []
for ep, path in files.items():
    df = pd.read_csv(path, sep="\t")
    df["epitope"] = ep
    dfs.append(df)

all_df = pd.concat(dfs, ignore_index=True)

# 2. تجهيز الـfeatures
# تحويل CDR3β إلى k-mer counts (مثلاً 3-mers)
vectorizer = CountVectorizer(analyzer="char", ngram_range=(3,3))
cdr3_features = vectorizer.fit_transform(all_df["CDR3b"].fillna(""))

# تحويل TRBV/TRBJ إلى one-hot
enc = OneHotEncoder()
vj_features = enc.fit_transform(all_df[["TRBV","TRBJ"]].fillna(""))

# دمج كل الـfeatures
import scipy.sparse as sp
X = sp.hstack([cdr3_features, vj_features])

# Labels = epitopes
y = all_df["epitope"]

# 3. Train/test split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 4. تدريب موديل تصنيف (Random Forest)
from sklearn.ensemble import RandomForestClassifier

clf = RandomForestClassifier(n_estimators=200, random_state=42, class_weight="balanced")
clf.fit(X_train, y_train)
y_pred = clf.predict(X_test)

from sklearn.metrics import classification_report
print(classification_report(y_test, y_pred, zero_division=0))
