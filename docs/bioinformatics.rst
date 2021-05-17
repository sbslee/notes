Bioinformatics
**************

Frequently used commands for Bioinformatics
===========================================

* To extract regions from a BED file:

    .. code-block:: console

        $ awk '{print $1":"$2"-"$3}' example.bed | sed 's/chr//g' > regions.list
