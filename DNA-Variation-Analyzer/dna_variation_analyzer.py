import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO


# ============================================================
# 1. GET INFORMATION FROM THE COMMAND LINE
# ============================================================

fasta_file = sys.argv[1]
reference_id = sys.argv[2]


# ============================================================
# 2. READ THE FASTA FILE
# ============================================================

records = list(SeqIO.parse(fasta_file, "fasta"))


# ============================================================
# 3. FIND THE REFERENCE SEQUENCE
# ============================================================

reference = None

for record in records:

    if record.id == reference_id:
        reference = record
        break


# Check if reference was found
if reference is None:
    print("Reference sequence not found.")
    sys.exit()


print()
print("DNA Variation Analyzer")
print("-----------------------")
print("Reference ID:", reference.id)
print("Reference length:", len(reference.seq))


# ============================================================
# 4. SELECT SAMPLES
# ============================================================

# Skip the first 5 sequences
samples = records[5:]


# Lists for storing results
all_mutations = []
compared_sample_ids = []
length_mismatch_samples = []


# ============================================================
# 5. COMPARE REFERENCE WITH ALL SAMPLES
# ============================================================

for sample in samples:

    if sample.id == reference.id:
        continue

    print()
    print("Sample ID:", sample.id)
    print("Sample length:", len(sample.seq))

    # Check sequence length
    if len(reference.seq) != len(sample.seq):

        print("Length is different. Skipping this sample.")

        length_mismatch_samples.append({
            "sample_id": sample.id,
            "sample_length": len(sample.seq),
            "reference_length": len(reference.seq)
        })

        continue

    # This sample was successfully compared
    compared_sample_ids.append(sample.id)

    mutations = []

    # Compare the two sequences
    for position in range(len(reference.seq)):

        reference_base = reference.seq[position]
        sample_base = sample.seq[position]

        if reference_base != sample_base:

            mutation = {
                "sample_id": sample.id,
                "position": position + 1,
                "reference": reference_base,
                "sample": sample_base
            }

            mutations.append(mutation)

    # Print a concise result
    if len(mutations) == 0:

        print("No mutations found.")

    else:

        print("Mutations found:", len(mutations))

    # Add mutations to the complete list
    all_mutations.extend(mutations)


print()
print("Analysis is complete")
print("TOTAL MUTATIONS:", len(all_mutations))


# ============================================================
# 6. CREATE MUTATION DATAFRAME
# ============================================================

df = pd.DataFrame(all_mutations)

print()
print("Mutation details:")
print(df)


# ============================================================
# 7. MUTATION SUMMARY
# ============================================================

mutation_summary = (
    df.groupby(
        ["position", "reference", "sample"]
    )
    .size()
    .reset_index(name="mutation_count")
)

mutation_summary["mutation_type"] = (
    mutation_summary["reference"]
    + " → "
    + mutation_summary["sample"]
)

mutation_summary = mutation_summary[
    [
        "position",
        "reference",
        "sample",
        "mutation_type",
        "mutation_count"
    ]
]

mutation_summary = mutation_summary.sort_values(
    "mutation_count",
    ascending=False
).reset_index(drop=True)


print()
print("Mutation Summary:")
print(mutation_summary)


# ============================================================
# 8. SAMPLE SUMMARY
# ============================================================

sample_summary = (
    df.groupby("sample_id")
    .size()
    .reset_index(name="mutation_count")
)


# Find valid samples that had zero mutations
zero_mutation_samples = [
    sample_id
    for sample_id in compared_sample_ids
    if sample_id not in sample_summary["sample_id"].values
]


# Add zero-mutation samples
zero_df = pd.DataFrame({
    "sample_id": zero_mutation_samples,
    "mutation_count": 0
})


sample_summary = pd.concat(
    [sample_summary, zero_df],
    ignore_index=True
)


sample_summary = sample_summary.sort_values(
    "mutation_count",
    ascending=False
).reset_index(drop=True)


print()
print("Sample Summary:")
print(sample_summary)


print()
print("Samples with zero mutations:",
      len(zero_mutation_samples))

print(zero_mutation_samples)


# ============================================================
# 9. POSITION SUMMARY
# ============================================================

position_summary = (
    df.groupby("position")["sample_id"]
    .nunique()
    .reset_index(name="sample_count")
)

position_summary = position_summary.sort_values(
    "sample_count",
    ascending=False
).reset_index(drop=True)


print()
print("Most Variable Positions:")
print(position_summary.head(20))


# ============================================================
# 10. MUTATION TYPE SUMMARY
# ============================================================

mutation_type_summary = (
    df.groupby(["reference", "sample"])
    .size()
    .reset_index(name="mutation_count")
)

mutation_type_summary["mutation_type"] = (
    mutation_type_summary["reference"]
    + " → "
    + mutation_type_summary["sample"]
)

mutation_type_summary = mutation_type_summary[
    [
        "mutation_type",
        "mutation_count"
    ]
].sort_values(
    "mutation_count",
    ascending=False
).reset_index(drop=True)


print()
print("Mutation Type Summary:")
print(mutation_type_summary)


# ============================================================
# 11. CREATE MASTER MUTATION TABLE
# ============================================================

master_table = mutation_summary.copy()


master_table = master_table.merge(
    position_summary,
    on="position",
    how="left"
)


master_table = master_table[
    [
        "position",
        "reference",
        "sample",
        "mutation_type",
        "mutation_count",
        "sample_count"
    ]
]


master_table = master_table.sort_values(
    "mutation_count",
    ascending=False
).reset_index(drop=True)


print()
print("Master Mutation Table:")
print(master_table)


# ============================================================
# 12. VALIDATION
# ============================================================

print()
print("--- Validation ---")


# Check total number of mutations
print(
    "Total mutations:",
    len(all_mutations)
)

print(
    "Total mutations from master table:",
    master_table["mutation_count"].sum()
)


# Check number of unique mutation types
print(
    "Original unique mutations:",
    df[
        ["position", "reference", "sample"]
    ].drop_duplicates().shape[0]
)

print(
    "Master table unique mutations:",
    len(master_table)
)


# Check sample counts
calculated_sample_count = (
    df.groupby("position")["sample_id"]
    .nunique()
    .reset_index(name="calculated_sample_count")
)


sample_count_check = master_table[
    ["position", "sample_count"]
].drop_duplicates().merge(
    calculated_sample_count,
    on="position",
    how="left"
)


sample_count_check["match"] = (
    sample_count_check["sample_count"]
    ==
    sample_count_check["calculated_sample_count"]
)


print()
print("All sample counts correct:",
      sample_count_check["match"].all())


# Check length mismatches
print()
print(
    "Length-mismatched samples:",
    len(length_mismatch_samples)
)


if length_mismatch_samples:

    length_mismatch_df = pd.DataFrame(
        length_mismatch_samples
    )

    print(length_mismatch_df)

else:

    length_mismatch_df = pd.DataFrame(
        columns=[
            "sample_id",
            "sample_length",
            "reference_length"
        ]
    )

    print("No length mismatches found.")


# ============================================================
# 13. EXPORT RESULTS TO EXCEL
# ============================================================

os.makedirs("results", exist_ok=True)
master_table.to_excel(
    "results/master_mutation_table.xlsx",
    index=False
)


mutation_summary.to_excel(
    "results/mutation_summary.xlsx",
    index=False
)


sample_summary.to_excel(
    "results/sample_summary.xlsx",
    index=False
)


position_summary.to_excel(
    "results/position_summary.xlsx",
    index=False
)


mutation_type_summary.to_excel(
    "results/mutation_type_summary.xlsx",
    index=False
)


df[
    [
        "sample_id",
        "position",
        "reference",
        "sample"
    ]
].to_excel(
    "results/mutation_details.xlsx",
    index=False
)


length_mismatch_df.to_excel(
    "results/length_mismatch_samples.xlsx",
    index=False
)


print()
print("Excel files created successfully.")


# ============================================================
# 14. CREATE FIGURES
# ============================================================

# -------- Mutation Type Distribution --------

os.makedirs("figures", exist_ok=True)
plt.figure()

plt.bar(
    mutation_type_summary["mutation_type"],
    mutation_type_summary["mutation_count"]
)

plt.xlabel("Mutation Type")
plt.ylabel("Number of Mutations")
plt.title("Distribution of Mutation Types")

plt.xticks(rotation=45)
plt.tight_layout()

plt.savefig(
    "figures/distribution_of_mutation.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# -------- Top 10 Samples --------

top_sample = sample_summary.head(10)


plt.figure()

plt.bar(
    top_sample["sample_id"],
    top_sample["mutation_count"]
)

plt.xlabel("Sample ID")
plt.ylabel("Number of Mutations")
plt.title("Top 10 Samples by Mutation Count")

plt.xticks(rotation=45)
plt.tight_layout()

plt.savefig(
    "figures/top_10_sample.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


print()
print("Figures created successfully.")