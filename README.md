pdbgrep
=======

Print atom records from Protein Data Bank (PDB) files that contain fields
that match a pattern.
 
Usage
-----

`pdbgrep [options] [file ...]`
 
Options
-------

    -a range        - atom number (ordinal)
    -b range        - temperature factor
    -c expression   - chain id
    -e expression   - element type
    -h              - heterogen atom
    -i expression   - insertion code
    -l expression   - alternate location indicator
    -n expression   - atom name
    -o range        - occupancy
    -p              - pass non-ATOM or HETATM records through to the output
    -q expression   - atomic charge
    -r range        - residue number
    -s expression   - segment id
    -t expression   - residue type
    -v              - left-shifted atom name (heavy atom or hydrogen)
    -x range        - x coordinate
    -y range        - y coordinate
    -z range        - z coordinate

*range*: a single numeric value, or `low:high`, or `low:` (no upper bound) or `:high` (no lower bound) 

*expression*: any regular expression that matches the field completely (leading and trailing spaces in fields are discarded)

The uppercase version of each option prints records that do NOT match the range or expression.

License
-------
MIT License

