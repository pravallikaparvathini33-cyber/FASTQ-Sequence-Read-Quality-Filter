# FASTQ Sequence Read Quality Filter

This project converts the original FASTQ quality-filtering Python program into an interactive Streamlit web application.

## Project features
- Uses Biopython `SeqIO` to parse FASTQ input.
- Calculates average per-base Phred quality score.
- Uses a user-adjustable minimum quality threshold (default: 20).
- Separates good-quality and low-quality reads.
- Calculates GC content for each read.
- Displays total, passed and rejected read counts.
- Displays read-level quality and GC-content results.
- Downloads passing reads and rejected reads as FASTQ files.
- Includes the eight example reads from the original project.

## Filtering logic
For each read:
1. Read the per-base Phred quality scores.
2. Calculate their average.
3. Calculate GC content as `(G + C) / sequence length * 100`.
4. If average quality is greater than or equal to the selected threshold, retain the read; otherwise reject it.

Default quality threshold: 20.

## Run locally
```bash
python -m pip install -r requirements.txt
streamlit run streamlit_app.py
```

## Deploy on Streamlit Community Cloud
1. Create a **public** GitHub repository.
2. Upload `streamlit_app.py`, `requirements.txt`, and `README.md`.
3. Visit https://share.streamlit.io/ and sign in with GitHub.
4. Create an app from the repository.
5. Set the main file path to `streamlit_app.py`.
6. Deploy and share the generated public `streamlit.app` URL.

## FASTQ sample-data note
The original sample text had two sequence/quality length mismatches:
- Read5 had a 48-base sequence but 50 quality characters.
- Read6 had a 49-base sequence but 48 quality characters.

Those two quality strings were adjusted in the built-in demo only so each FASTQ record is valid and Biopython can parse it. The sequence strings and intended quality categories are preserved. For a real submission, you can upload your own valid FASTQ file instead.
