import streamlit as st
from Bio import SeqIO
from io import StringIO
import pandas as pd

# ==========================================
# FASTQ SEQUENCE READ QUALITY FILTER
# Biopython + Streamlit Web Application
# ==========================================

st.set_page_config(
    page_title="FASTQ Sequence Read Quality Filter",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 FASTQ Sequence Read Quality Filter")
st.write(
    "A Biopython web application that filters FASTQ reads using "
    "average Phred quality score and reports GC content."
)

st.sidebar.header("Filtering Settings")
quality_threshold = st.sidebar.slider(
    "Minimum average quality score",
    min_value=0,
    max_value=41,
    value=20,
    help="A read is retained when its average Phred quality is at least this value."
)

# The sample reads from the supplied Python project.
# Quality strings are kept the same length as their sequences so the
# example is valid FASTQ and can be parsed by Biopython.
BUILT_IN_FASTQ = """@Read1
ATGCGTACGTAGCTAGCTAGATCGATCGTAGCTAGCTAGCTAGCTAGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@Read2
ATGCTAGCTAGCTAACGTACGTAGGCTACGATCGTACGTAGCTAGCTA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@Read3
ATGCGCGATATATCGCGCTATGCGCGATATCGATCGCGTATATCGCGC
+
555555555555555555555555555555555555555555555555
@Read4
GCTAGCTAGCATCGATCGATGCTAGCTAGCATCGATCGATGCTAGCTA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@Read5
ATATATATATATATATATATATATATATATATATATATATATATATAT
+
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
@Read6
CGTACGTAGCTAGCTAGGCTACGTAGCTAGCTACGTAGCTAGCTAGCTA
+
5555555555555555555555555555555555555555555555555
@Read7
ATGCGTACGTAGCTAGCTAGATCGATCGTAGCTAGCTAGCTAGCTAGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@Read8
GCGCGCGCGCATATATATATGCGCGCGCATATATATATGCGCGCGCAT
+
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
"""

st.subheader("1. Select FASTQ input")
input_choice = st.radio(
    "Choose how to provide your FASTQ data:",
    ["Use the project's built-in example reads", "Upload my own FASTQ file"],
    horizontal=True
)

fastq_text = None

if input_choice == "Use the project's built-in example reads":
    fastq_text = BUILT_IN_FASTQ
    st.info("Using the 8 example reads from your original Python project.")
else:
    uploaded_file = st.file_uploader(
        "Upload a FASTQ file",
        type=["fastq", "fq"],
        help="Upload a standard uncompressed FASTQ file."
    )
    if uploaded_file is not None:
        try:
            fastq_text = uploaded_file.getvalue().decode("utf-8")
        except UnicodeDecodeError:
            st.error("Could not read this file as text. Please upload a standard uncompressed FASTQ file.")

if fastq_text is not None:
    try:
        # ------------------------------------------
        # READ FASTQ DATA USING BIOPYTHON
        # ------------------------------------------
        reads = list(SeqIO.parse(StringIO(fastq_text), "fastq"))

        if not reads:
            st.error("No valid FASTQ reads were found in this input.")
            st.stop()

        good_reads = []
        bad_reads = []
        result_rows = []

        # ------------------------------------------
        # QUALITY SCORE + GC CONTENT CALCULATION
        # ------------------------------------------
        for record in reads:
            sequence = str(record.seq).upper()
            quality_scores = record.letter_annotations["phred_quality"]

            avg_quality = (
                sum(quality_scores) / len(quality_scores)
                if quality_scores else 0.0
            )

            gc_count = sequence.count("G") + sequence.count("C")
            gc_content = (
                (gc_count / len(sequence)) * 100
                if sequence else 0.0
            )

            if avg_quality >= quality_threshold:
                good_reads.append(record)
                status = "Passed"
            else:
                bad_reads.append(record)
                status = "Rejected"

            result_rows.append({
                "Read ID": record.id,
                "Sequence Length (bp)": len(sequence),
                "Average Quality": round(avg_quality, 2),
                "GC Content (%)": round(gc_content, 2),
                "Result": status
            })

        # ------------------------------------------
        # DISPLAY FILTERING RESULTS
        # ------------------------------------------
        st.subheader("2. Filtering results")
        total_reads = len(reads)
        good_count = len(good_reads)
        bad_count = len(bad_reads)
        pass_percentage = (good_count / total_reads) * 100

        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Total reads", total_reads)
        col2.metric("Good quality reads", good_count)
        col3.metric("Low quality reads", bad_count)
        col4.metric("Reads passed", f"{pass_percentage:.1f}%")

        st.progress(pass_percentage / 100)
        st.caption(
            f"Filtering rule: average Phred quality score ≥ {quality_threshold}"
        )

        # ------------------------------------------
        # READ-BY-READ DETAILS
        # ------------------------------------------
        st.subheader("3. Read details")
        results_df = pd.DataFrame(result_rows)
        st.dataframe(results_df, use_container_width=True, hide_index=True)

        st.subheader("4. Good quality reads")
        if good_reads:
            for record in good_reads:
                sequence = str(record.seq).upper()
                scores = record.letter_annotations["phred_quality"]
                avg_quality = sum(scores) / len(scores) if scores else 0
                gc_count = sequence.count("G") + sequence.count("C")
                gc_content = (gc_count / len(sequence)) * 100 if sequence else 0

                with st.expander(f"{record.id} — Average quality: {avg_quality:.2f}"):
                    st.write("**Sequence:**")
                    st.code(sequence)
                    st.write(f"**Average Quality:** {avg_quality:.2f}")
                    st.write(f"**GC Content:** {gc_content:.2f}%")
        else:
            st.warning("No reads passed the selected quality threshold.")

        # ------------------------------------------
        # FILTERED FASTQ DOWNLOAD
        # ------------------------------------------
        st.subheader("5. Download filtered FASTQ")
        if good_reads:
            output = StringIO()
            SeqIO.write(good_reads, output, "fastq")
            filtered_fastq = output.getvalue()

            st.download_button(
                "⬇️ Download good quality reads (FASTQ)",
                data=filtered_fastq,
                file_name="filtered_good_quality_reads.fastq",
                mime="text/plain"
            )
        else:
            st.info("There are no passing reads to download.")

        # Optional rejected-read download
        if bad_reads:
            rejected_output = StringIO()
            SeqIO.write(bad_reads, rejected_output, "fastq")
            st.download_button(
                "Download rejected reads (FASTQ)",
                data=rejected_output.getvalue(),
                file_name="rejected_low_quality_reads.fastq",
                mime="text/plain"
            )

        st.success("FASTQ quality filtering completed successfully.")

    except Exception as error:
        st.error(f"Unable to process this FASTQ input: {error}")
        st.info(
            "Please check that the file is valid FASTQ. Each read must have "
            "a header, sequence, plus line and quality string of matching length."
        )
else:
    st.info("Choose the built-in example or upload a FASTQ file to begin.")

st.divider()
st.caption(
    "Project: FASTQ Sequence Read Quality Filter | "
    "Python • Biopython • Streamlit"
)
