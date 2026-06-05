
import pandas as pd
import json

# 1. قراءة الملف الأصلي
# تأكد أن اسم الملف صح أو غيره لاسم ملفك
input_file = 'VDJdb_KTWGQYWQV_TRA_TRB.tsv'
df = pd.read_csv(input_file, sep='\t')

# 2. وظيفة لتنظيف أسماء الجينات (إزالة *01 وما بعدها)
def clean_gene(gene_name):
    if pd.isna(gene_name):
        return ""
    return str(gene_name).split('*')[0].split('/')[0] # بياخد أول جزء قبل النجمة أو السلاش

# 3. وظيفة لاستخراج الـ Subject ID من عمود Meta
def extract_subject(meta_str):
    try:
        # تحويل النص لقاموس (Dictionary)
        meta_dict = json.loads(meta_str.replace("'", '"'))
        # البحث عن المعرف، لو مش موجود نستخدم الكوهورت، لو مش موجود نكتب "Unknown"
        subject = meta_dict.get('subject.id') or meta_dict.get('subject.cohort') or "Sub1"
        return subject
    except:
        return "Sub1"

# 4. فصل وتحضير بيانات الـ Beta
beta = df[df['Gene'] == 'TRB'].copy()
beta['subject'] = beta['Meta'].apply(extract_subject)
beta['TRBV'] = beta['V'].apply(clean_gene)
beta['TRBJ'] = beta['J'].apply(clean_gene)
beta = beta[['complex.id', 'CDR3', 'TRBV', 'TRBJ', 'subject']]
beta.columns = ['complex.id', 'CDR3b', 'TRBV', 'TRBJ', 'subject']

# 5. فصل وتحضير بيانات الـ Alpha
alpha = df[df['Gene'] == 'TRA'].copy()
alpha['TRAV'] = alpha['V'].apply(clean_gene)
alpha['TRAJ'] = alpha['J'].apply(clean_gene)
alpha = alpha[['complex.id', 'CDR3', 'TRAV', 'TRAJ']]
alpha.columns = ['complex.id', 'CDR3a', 'TRAV', 'TRAJ']

# 6. دمج الجدولين في سطر واحد بناءً على complex.id
final_df = pd.merge(beta, alpha, on='complex.id', how='left')

# 7. إضافة عمود الـ count (مطلوب لـ GLIPH2)
final_df['count'] = 1

# 8. ترتيب الأعمدة النهائي كما في الصورة المقبولة
final_columns = ['CDR3b', 'TRBV', 'TRBJ', 'CDR3a', 'TRAV', 'TRAJ', 'subject', 'count']
final_df = final_df[final_columns]

# 9. حفظ الملف بصيغة TAB (أهم خطوة للقبول)
output_filename = 'VDJdb_KTWGQYWQV_TRA_TRB_Input.csv'
final_df.to_csv(output_filename, sep='\t', index=False)

print(f"تم بنجاح! الملف جاهز للرفع على GLIPH2 باسم: {output_filename}")