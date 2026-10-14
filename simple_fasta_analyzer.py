# simple_fasta_analyzer.py
# A beginner-friendly FASTA analyzer using only built-in Python (no imports).


def read_fasta(filepath):
    """Read a FASTA file and return a list of (header, sequence) pairs."""
    sequences = []
    header = None
    seq_parts = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue  # skip empty lines

            if line.startswith(">"):
                # A new header: save the previous sequence first
                if header is not None:
                    sequences.append((header, "".join(seq_parts)))
                header = line[1:].strip()
                seq_parts = []
            else:
                if header is None:
                    raise ValueError("File must start with a '>' header line.")
                seq_parts.append(line.replace(" ", "").upper())

    # Don't forget the last sequence in the file
    if header is not None:
        sequences.append((header, "".join(seq_parts)))

    if not sequences:
        raise ValueError("No sequences found in the file.")
    return sequences


def count_bases(seq):
    """Count A, C, G, T, N. Any other symbol is counted as 'other'."""
    counts = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0, "other": 0}
    for base in seq:
        if base in "ACGTN":
            counts[base] += 1
        else:
            counts["other"] += 1
    return counts


def gc_content(counts):
    """GC% = (G + C) / (A + C + G + T) * 100. Unknown bases are not included."""
    known = counts["A"] + counts["C"] + counts["G"] + counts["T"]
    if known == 0:
        return 0.0
    return (counts["G"] + counts["C"]) / known * 100


def analyze(sequences):
    """Build one result (a dictionary) for each sequence."""
    results = []
    for header, seq in sequences:
        if not seq:
            print("Warning: '{}' has no sequence, skipped.".format(header))
            continue
        counts = count_bases(seq)
        results.append({
            "header": header,
            "length": len(seq),
            "gc": gc_content(counts),
            "counts": counts,
        })
    return results


def print_table(results):
    """Print the results as a neat table."""
    row = "{:<40} {:>8} {:>8} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6}"
    title = row.format("Header", "Length", "GC%", "A", "C", "G", "T", "N", "Other")
    print(title)
    print("-" * len(title))  # separator always matches the table width

    for r in results:
        c = r["counts"]
        name = r["header"][:40]  # cut long headers so columns stay aligned
        print(row.format(name, r["length"], "{:.2f}".format(r["gc"]),
                         c["A"], c["C"], c["G"], c["T"], c["N"], c["other"]))


def save_csv(results, csv_path):
    """Save the results to a CSV file."""
    with open(csv_path, "w") as out:
        out.write("header,length,gc_percent,A,C,G,T,N,other\n")
        for r in results:
            c = r["counts"]
            name = r["header"].replace(",", ";")  # commas would break the CSV
            out.write("{},{},{:.2f},{},{},{},{},{},{}\n".format(
                name, r["length"], r["gc"],
                c["A"], c["C"], c["G"], c["T"], c["N"], c["other"]))


def make_output_path(fasta_path):
    """'data/sample.fasta' -> 'data/sample_summary.csv'"""
    dot = fasta_path.rfind(".")
    slash = max(fasta_path.rfind("/"), fasta_path.rfind("\\"))
    if dot > slash:                 # the dot is in the file name, not a folder
        fasta_path = fasta_path[:dot]
    return fasta_path + "_summary.csv"


def main():
    fasta_path = input("Enter path to FASTA file: ").strip()
    if not fasta_path:
        print("No file path provided. Exiting.")
        return

    try:
        sequences = read_fasta(fasta_path)
        results = analyze(sequences)
        if not results:
            print("No valid sequences to analyze.")
            return

        print("\nSummary of sequences:")
        print_table(results)

        csv_path = make_output_path(fasta_path)
        save_csv(results, csv_path)
        print("\nSummary saved to:", csv_path)

    except FileNotFoundError:
        print("Error: file not found:", fasta_path)
    except ValueError as e:
        print("Error:", e)


if __name__ == "__main__":
    main()
