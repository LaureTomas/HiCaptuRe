---
editor_options: 
  markdown: 
    wrap: 72
---

## CHANGES IN VERSION 0.10.0

New class HiCaptuRe

New function load_genome and digest_genome

New message when loading the package

Renamed columns gene\_\* by bait\_\*

Renamed in column names \_I and \_II by \_1 and \_2

Renamed P (promoter) by B (bait)

Include PAR file for Homo Sapiens GRCh38

### Class HiCapture

-   Inherits methods from GenomicInteractions
-   New slots
    -   parameters: include parameters used in different functions
        applied to the object
    -   ByBait: list containing bait centric statistics
    -   ByRegion: list containing region centric statistics

### load_interactions

-   Accepts bedpe files
-   Returns HiCaptuRe object
-   seqinfo updated based on digest genome
-   Interactions are sorted increasingly by genomic regions of both ends
-   Added several errors if genome, or digest not correct
-   Show progress bar with progressr::handlers(global = T) and extra
    messages with progressr::handlers("progress")
-   New argument washU_seqname
-   Extra columns
    -   ID1 and ID2: identifier of restriction fragment based on digest
        genome
    -   distance

### annotate_interactions

-   Remove UCE annotation
-   Reannote the column int based on the new annotation of fragments
-   Internal function annotate_POEuce renamed to annotate_BOE
-   Reorder interactions based on ID, in case some fragments are now
    annotated as baits

### export_interactions

-   Export parameters file
-   When exporting HiCaptuRe object with more than 1 CS column it will
    generate several files

### distance_summary

-   Column ibed renamed to HiCaptuRe

### interactionsByBaits

-   New argument sep to indicate separator in bait names
-   New slot ByBaits containing Bait centric information

### interactionsByRegions

-   Removed nomenclature \_II or \_I... to indicate the presence of an
    overlaping region in any end
-   Add several new columns regarding overlaping peaks:
    -   region: T/F overlaping region
    -   Nregion: number of overlaping regions
    -   regionID: ID of overlaping regions
    -   regionCov: coverage of the overlapoing regions on that fragment
-   New slot ByRegions containing Region centric information
