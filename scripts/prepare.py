#!/usr/bin/env python3
#usage: python3 prepare.py -i metadata_small_final.csv

# imports:
import os
import argparse
import pandas as pd


def read_file(file):
    """ Read the metadata file containing path files and return a pandas dataframe"""
    metadata = pd.read_csv(file, sep='\t', dtype=str)
    return metadata


def check_mapdirs(mapdirs):
    """ Check if the mapdir exists or not"""
    for mapdir in mapdirs:
        if not os.path.exists(mapdir):
            raise Exception("WARNING following mapdir is missing %s" % mapdir)


def unique(liste):
    return sorted(set(liste))


def get_path_files(df):
    """ Group the samples by "sample".
    Return a dataframe with one line per individual.
    """
    df = df.groupby('sample').agg(lambda x: unique(list(x))).reset_index()
    return df


def get_species_short(species):
    animal_keys = ["bos_taurus", "sus_scrofa"]
    animal_values = ["cow", "pig"]
    succint_name = dict(zip(animal_keys, animal_values))
    if species in succint_name:
        return succint_name[species]
    else:
        return species


def write_sample_file(metadata, output_file):
    """Read the metadate file and return a pandas object"""
    sample_ids = []
    species_names = []
    for index, row in metadata.iterrows():
        short_name = get_species_short(row['species'])
        sample = "-".join([short_name, row['tissue'],
                           row['animal_name']])
        species = short_name
        sample_ids.append(sample)
        species_names.append(species)
    metadata['sample'] = sample_ids
    metadata['species'] = species_names
    samples = metadata[['sample', 'sample_unit', 'species', 'sex', 'mapdir']]
    samples.to_csv(output_file, sep="\t", index=False)


def main():

    # Read the metadata file:
    metadata = read_file(args.input_file)

    # Check the mapdir of each sample exists or not:
    check_mapdirs(metadata["mapdir"])

    # Write a new output file 'samples.tsv' with identical key per individual:
    write_sample_file(metadata, args.output_file)


def parse_arguments():
    parser = argparse.ArgumentParser(description='Metadata file')
    parser.add_argument('-i', '--input_file',
                        required=True, help='Metadata file')
    parser.add_argument('-o', '--output_file',
                        required=True, help='Metadata file')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    main()
