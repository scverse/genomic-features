from genomic_features.ensembl.ensembldb import EnsemblDB


def listColumns(ensdb: EnsemblDB, table=None) -> list:
    """
    List all columns available in the genomic features table.
    """

    if table is None:
        table = ensdb.db.list_tables()  # list of table names
    else:
        table = [table]  # list of tables names (only one)
    columns = []
    for t in table:
        columns += ensdb.db.table(t).columns
    return columns


def cleanColumns(ensdb: EnsemblDB, columns) -> list:
    """
    Clean a list of columns to make sure they are valid.
    """
    valid_columns = listColumns(ensdb)
    columns = [c for c in columns if c in valid_columns]
    unvalid_columns = [c for c in columns if c not in valid_columns]

    if len(unvalid_columns) > 0:
        print("The following columns are not valid and will be ignored:")
        print(unvalid_columns)

    if len(columns) == 0:
        raise ValueError("No valid columns were found.")
    return columns


def tablesForColumns(ensdb: EnsemblDB, columns) -> list:
    """
    Return a list of tables that contain the specified columns.
    """
    columns = cleanColumns(ensdb, columns)
    tables = []
    for t in ensdb.db.list_tables():
        if set(columns).issubset(ensdb.db.table(t).columns):
            tables.append(t)
    return tables
