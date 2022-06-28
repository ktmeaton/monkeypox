import pandas as pd
import argparse
import numpy as np

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="remove time info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, required=True, help="input data")
    parser.add_argument('--strain-id', type=str, required=True, help="field to use as strain id")
    parser.add_argument('--output', type=str, required=True, help="output metadata")
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, sep='\t')
    if 'strain' in metadata.columns:
        metadata.rename(columns={'strain': 'strain_original'}, inplace=True)

    # copy column, retaining original
    # ie keep "accession" but also include "strain" with data from previous "accession" column
    # insert this as the first column in the dataframe
    metadata.insert(0, "strain", metadata[args.strain_id])
    # metadata["strain"] = metadata[args.strain_id]

    # Remove the year column, because it will break augur filter
    if "year" in metadata.columns:
        new_dates = []
        # Iterate through the `date`` and `year`` fields
        for s_date, s_year in zip(metadata["date"], metadata["year"]):

            # If date is null, we use the year
            if pd.isna(s_date) and not pd.isna(s_year):
                new_dates.append("{}-XX-XX".format(int(s_year)))

            # if date is not null, use it
            elif not pd.isna(s_date):
                new_dates.append(s_date)

            # Otherwise, use none
            else:
                new_dates.append(None)

        metadata["date"] = new_dates
        metadata.drop(columns=["year"], inplace=True)
        
    metadata.to_csv(args.output, sep='\t', index=False)
