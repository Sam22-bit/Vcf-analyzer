# ====== ABSOLUTE FIRST COMMANDS ======
import streamlit as st
st.set_page_config(
    page_title="VCF Analyzer",
    layout="wide",
    page_icon="🧬",
    initial_sidebar_state="expanded"
)

background_url = 'https://img.freepik.com/premium-vector/dna-strand-background-genetic-engineering-laboratory-research_7117149.jpg'
page_bg_img = f'''
<style>
[data-testid="stAppViewContainer"] > .main {{
    background-image: url("{background_url}")  !important;
    background-size: cover;
    background-position: center;
    background-repeat: no-repeat;
    background-attachment: fixed;
}}

[data-testid="stHeader"] {{
    background: rgba(0,0,0,0);
}}

[data-testid="stToolbar"] {{
    right: 2rem;
}}
</style>
'''

st.markdown(page_bg_img, unsafe_allow_html=True)

# ====== THEME CONFIGURATION ======
st.markdown("""<style>
/* Main theme colors - light theme by default */
:root {
    --primary-color: #f3cbe5;
    --background-color: #fbfbfb;
    --secondary-background-color: #e8ecfd;
    --text-color: #000000;
}

/* Force light theme */
[data-testid="stTheme"] {
    color-scheme: light;
}

/* Apply theme to components */
.stApp {
    background-color: var(--background-color);
}
</style>
""", unsafe_allow_html=True)

# ====== IMPORTS ======
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tempfile
import os
from collections import Counter
from cyvcf2 import VCF

# --- Enhanced VCF Parser ---
def parse_vcf(vcf_path):
    try:
        vcf = VCF(vcf_path)
        records = []
        for variant in vcf:
            # Get NGS-specific metrics if available
            alt_count = 0
            total_depth = 0
            if 'AD' in variant.FORMAT:
                sample = variant.genotypes[0]  # First sample
                ref_count, alt_count = sample['AD'][0], sample['AD'][1]
                total_depth = ref_count + alt_count
                vaf = alt_count / total_depth if total_depth > 0 else 0
            else:
                vaf = variant.INFO.get('AF', [0])[0] if 'AF' in variant.INFO else 0

            records.append({
                'CHROM': variant.CHROM,
                'POS': variant.POS,
                'REF': variant.REF,
                'ALT': ','.join(variant.ALT),
                'TYPE': variant.var_type,
                'AF': variant.INFO.get('AF'),
                'VAF': vaf,
                'DP': total_depth,
                'Mutation': f"{variant.REF}→{variant.ALT[0]}",
                'QUAL': variant.QUAL
            })
        return pd.DataFrame(records)
    except Exception as e:
        st.error(f"VCF parsing error: {str(e)}")
        return pd.DataFrame()

# --- Analysis Functions ---
def calculate_mutation_frequencies(df):
    if df.empty:
        return pd.DataFrame()

    mutation_counts = Counter(df['Mutation'])
    total_mutations = sum(mutation_counts.values())

    freq_df = pd.DataFrame.from_dict(mutation_counts, orient='index', columns=['Count'])
    freq_df['Frequency (%)'] = (freq_df['Count'] / total_mutations) * 100
    return freq_df.sort_values('Count', ascending=False)

def ngs_quality_metrics(df):
    metrics = {}
    if not df.empty:
        metrics['Mean Depth'] = np.mean(df['DP'].dropna())
        metrics['Median VAF'] = np.median(df['VAF'].dropna()) * 100
        metrics['Q30 Variants'] = len(df[df['QUAL'] >= 30])
    return metrics

def mutation_spectrum(df):
    ti, tv = 0, 0
    transitions = [('A','G'), ('G','A'), ('C','T'), ('T','C')]
    transversions = [
        ('A','C'), ('A','T'), ('C','A'), ('C','G'),
        ('G','C'), ('G','T'), ('T','A'), ('T','G')
    ]

    if df.empty or 'snp' not in df['TYPE'].values:
        return {'Transitions': 0, 'Transversions': 0, 'Ti/Tv Ratio': 'NA'}
    
    for _, row in df[df['TYPE'] == 'snp'].iterrows():
        ref, alt = row['REF'], row['ALT'].split(',')[0]
        if (ref, alt) in transitions:
            ti += 1
        elif (ref, alt) in transversions:
            tv += 1
    return {
        'Transitions': ti, 
        'Transversions': tv, 
        'Ti/Tv Ratio': round(ti/tv, 2) if tv > 0 else 'NA'
    }

# --- Streamlit UI ---
# Sidebar for navigation
page = st.sidebar.radio("Explore", ("🏠 Home", "ℹ️ About", "🧬 VCF Analyzer"))

if page == "🏠 Home":
    st.markdown("""
    <style>
    .welcome-container {
        display: flex;
        flex-direction: column;
        align-items: center;
        justify-content: center;
        text-align: center;
    }
    </style>
    
    <div class="welcome-container">
        <h1 style='text-align: center;'>🎉 Welcome to the VCF Analyzer! 🧬</h1>
        <img src='https://media.giphy.com/media/ASd0Ukj0y3qMM/giphy.gif' width='250'>
        <h2>👋 Hello!</h2>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    **🔍 What is this?**

    - A simple yet **powerful web app** for analyzing Variant Call Format (VCF) files.
    - Built with ❤️ using **Streamlit**, **cyvcf2**, and **Matplotlib**.

    **🛠️ How to Use:**

    - Upload your `.vcf` or `.vcf.gz` file on the 'VCF Analyzer' page.
    - Instantly get insights: variant stats, quality metrics, and mutation spectra.

    **🚀 Use Cases:**

    - Quick exploration of VCF data.
    - Understanding NGS quality and variant distribution.
    - Visualizing allele frequencies and mutation patterns.

    **⚙️ What's Under the Hood:**

    - **VCF Parsing:** cyvcf2 for fast variant extraction.
    - **Visualization:** Matplotlib + Streamlit native charts.
    - **Data Handling:** Pandas + Numpy magic.

    **✨ Highlight Features:**

    - 📊 Variant Type Pie Chart with adaptive labeling.
    - 🗺️ Chromosome-wise distribution bar chart.
    - 🔬 Mutation Spectrum analysis with Ti/Tv ratio.
    - 📈 Allele Frequency histograms.

    **🔮 Future Goals:**

    - Add annotation with ClinVar / dbSNP.
    - Support for multi-sample VCFs.
    - Exportable detailed reports.
    """)

elif page == "ℹ️ About":
    st.header("👤 Our Team")
    col1, = st.columns(1)
    with col1:
        st.subheader("Author")
        st.image("https://media.licdn.com/dms/image/v2/D4E03AQF0-nTDDZexCg/profile-displayphoto-shrink_400_400/B4EZa1t2.nGcAg-/0/1746805421056?e=1752105600&v=beta&t=xdSUVzGH_d13G28MPE_hNZ-wo8fN_s9hdfDmrcnJbkw", width=150)
        st.markdown("Samiksha Nitin Pasalkar")
        st.markdown("I am Samiksha Nitin Pasalkar, currently pursuing an M.Sc. in Bioinformatics at DES Pune University. This application is a result of both my academic curriculum and my personal interest in exploring data-driven solutions for biological research.The VCF Analyzer was created to address the challenges faced by researchers and students in interpreting and visualizing Variant Call Format (VCF) files, which are widely used in genomics to store gene sequence variations. Manually parsing and understanding these files can be complex and time-consuming. This app simplifies the process by providing easy-to-understand visualizations, including variant type distributions and mutation spectra, along with a preview of the uploaded data. It is intended as an educational and practical tool for anyone working with genetic variation data, enhancing both learning and research efficiency.")
        st.markdown("Linkedin:[https://www.linkedin.com/in/samiksha-pasalkar-b2bb39320/](https://www.linkedin.com/in/samiksha-pasalkar-b2bb39320/)")

        st.subheader("Mentor")
        st.image("https://media.licdn.com/dms/image/v2/D5603AQF9gsU7YBjWVg/profile-displayphoto-shrink_400_400/B56ZZI.WrdH0Ag-/0/1744981029051?e=1752105600&v=beta&t=F4QBDSEgjUvnBS00xPkKqPTLI0jQaMpYefaOzARY1Yg", width=150)
        st.markdown("Dr.Kushagra Kashayp")
        st.markdown("Assistant Professor (Bioinformatics), Department of Life Sciences, School of Science and Mathematics, DES Pune University")
        st.markdown("Special thanks to Dr. Kushagra Kashyap Professor at DES Pune University.This project is developed under the guidance of Dr. Kushagra Kashyap, my professor at the university, whose teaching in bioinformatics has been instrumental in shaping my understanding and passion for the subject.")
        st.markdown("Linkedin:https://www.linkedin.com/in/dr-kushagra-kashyap-b230a3bb")

elif page == "🧬 VCF Analyzer":
    st.title("🧬 VCF Analyzer")
    st.markdown("""
    Analyze VCF files with:
    - **Basic variant statistics**
    - **NGS quality metrics**
    - **Mutation frequency analysis**
    - **Variant allele frequencies**
    - **Mutation spectrum**
    """)

    uploaded_vcf = st.file_uploader("Upload VCF File", type=["vcf", "vcf.gz"])

    if uploaded_vcf:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".vcf") as tmp:
            tmp.write(uploaded_vcf.getvalue())
            tmp_path = tmp.name
        
        try:
            with st.spinner("Analyzing VCF..."):
                df = parse_vcf(tmp_path)
                
                if df.empty:
                    st.error("No variants detected!")
                    st.stop()
                
                # --- Basic Statistics ---
                st.header("🔎 Basic Statistics")
                col1, col2 = st.columns(2)
                with col1:
                    st.metric("Total Variants", len(df))
                    st.metric("Chromosomes", df['CHROM'].nunique())
                with col2:
                    variant_counts = df['TYPE'].value_counts()
                    st.metric("SNPs", variant_counts.get('snp', 0))
                    st.metric("INDELs", variant_counts.get('indel', 0))

                # --- Variant Types ---
                st.header("🧬 Variant Type Distribution", divider='rainbow') 

                # Prepare data - count types and handle small percentages
                type_counts = df['TYPE'].value_counts()
                total = type_counts.sum()

                # Combine small categories (<5%) into "Other" if needed
                if (type_counts / total * 100).lt(5).any():
                    threshold = total * 0.05
                    main_types = type_counts[type_counts >= threshold]
                    other_count = type_counts[type_counts < threshold].sum()
                    if other_count > 0:
                        type_counts = pd.concat([main_types, pd.Series({'Other': other_count})])

                # Create figure with constrained layout
                fig1, ax1 = plt.subplots(figsize=(6, 6), constrained_layout=True)

                # Custom autopct function to handle small percentages
                def autopct_format(pct):
                    return f'{pct:.1f}%' if pct >= 5 else ''

                # Plot with improved parameters
                wedges, texts, autotexts = ax1.pie(
                    type_counts,
                    labels=type_counts.index,
                    autopct=autopct_format,
                    startangle=90,
                    pctdistance=0.8,
                    labeldistance=1.1,
                    colors=plt.cm.Pastel1.colors,
                    wedgeprops={'linewidth': 1, 'edgecolor': 'white'},
                    textprops={'fontsize': 10}
                )

                # Improve label positioning
                plt.setp(autotexts, size=10, weight='bold')
                plt.setp(texts, size=10) 

                # Equal aspect ratio ensures pie is drawn as circle
                ax1.axis('equal')

                # Add legend for small slices if we created "Other"
                if 'Other' in type_counts:
                    ax1.legend(
                        wedges,
                        [f'{l} ({s/total*100:.1f}%)' for l, s in zip(type_counts.index, type_counts)],
                        title="Variant Types",
                        loc="center left",
                        bbox_to_anchor=(1, 0, 0.5, 1)
                    )

                # Display in Streamlit
                st.pyplot(fig1)

                # --- Chromosome Distribution ---
                st.header("🗺️ Variants per Chromosome")
                st.bar_chart(df['CHROM'].value_counts().sort_index())

                # --- Mutation Frequency Analysis ---
                st.header("🔢 Mutation Frequencies")
                freq_df = calculate_mutation_frequencies(df)
                
                col1, col2 = st.columns(2)
                with col1:
                    st.dataframe(freq_df.head(20))
                with col2:
                    fig2, ax2 = plt.subplots(figsize=(8, 4))
                    freq_df.head(10).sort_values('Count').plot(
                        kind='barh', y='Count', ax=ax2, color='teal'
                    )
                    ax2.set_xlabel("Occurrence Count")
                    st.pyplot(fig2)

                # --- NGS Quality Metrics ---
                st.header("📊 NGS Quality Metrics")
                metrics = ngs_quality_metrics(df)
                
                cols = st.columns(3)
                cols[0].metric("Mean Depth", f"{metrics.get('Mean Depth', 0):.1f}x")
                cols[1].metric("Median VAF", f"{metrics.get('Median VAF', 0):.1f}%")
                cols[2].metric("Q30 Variants", metrics.get('Q30 Variants', 0))

                # --- Allele Frequency ---
                st.header("📈 Allele Frequency Distribution")
                if df['AF'].notnull().sum() > 0:
                    fig3, ax3 = plt.subplots(figsize=(8, 4))
                    ax3.hist(df['AF'].dropna(), bins=20, color='skyblue', edgecolor='black')
                    ax3.set_xlabel('Allele Frequency')
                    ax3.set_ylabel('Count')
                    st.pyplot(fig3)
                else:
                    st.info("No Allele Frequency (AF) info found")

                # --- Mutation Spectrum ---
                st.header("🔬 Mutation Spectrum")
                spectrum = mutation_spectrum(df)
                cols = st.columns(3)
                cols[0].metric("Transitions", spectrum['Transitions'])
                cols[1].metric("Transversions", spectrum['Transversions'])
                cols[2].metric("Ti/Tv Ratio", spectrum['Ti/Tv Ratio'])
                
                fig4, ax4 = plt.subplots(figsize=(8, 4))
                ax4.bar(['Transitions', 'Transversions'], 
                       [spectrum['Transitions'], spectrum['Transversions']], 
                       color=['#4CAF50', '#FF5722'])
                ax4.set_ylabel('Count')
                st.pyplot(fig4)

                # --- Raw Data ---
                with st.expander("📄 Show Raw Data"):
                    st.dataframe(df)

        finally:
            os.unlink(tmp_path)
    else:
        st.info("👈 Please upload a VCF file to start analysis.")
