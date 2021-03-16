## Dependencies

KMC3, stxxl, sdsl-lite.

## Compiling

```
# Download dependencies
git submodule init
git submodule update

# Build KMC
cd KMC
make
cd ..

# Build sdsl-lite
cd sdsl-lite
sh install.sh
cd ..

# Build stxxl
cd stxxl
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DUSE_GNU_PARALLEL=ON -DCMAKE_INSTALL_PREFIX=./install
make
make install
cd ../..

# Build bufBOSS
make all

```

## Usage

There are three programs: `bufboss_build`, `bufboss_update` and `bufboss_query`. They will be compiled to the directory `./bin`. 

### Example

We recommend building the index out of a KMC database. For example:

```
MY_INPUT=data/reads.fna
K=31
./KMC/bin/kmc -v -k$K -m1 -ci1 -cs1 -fm $MY_INPUT temp/kmc_db temp
./bin/bufboss_build -o my_index -t temp --KMC temp/kmc_db
```

To update the index by adding the k-mers in a file and deleting the k-mers in another file, run the following:

```
MY_ADDITIONS=data/additions.fna 
MY_DELETIONS=data/deletions.fna 
./bin/bufboss_update -i my_index/ -o my_index --add $MY_ADDITIONS --del $MY_DELETIONS --buffer-fraction 0.1
```

This overwrites the previous index. The value of --buffer-fraction determines how often the buffer is flushed. If the buffer fraction is t, then the buffer is flushed when it has 10% of the edges compared to the static part of the data structure. To query sequences agaisnt the index, run the following:

```
MY_QUERIES=data/queries.fna
./bin/bufboss_query -i my_index -q $MY_QUERIES
```

This prints bits 1 and 0 in ASCII such that for every input sequence R, we print one line L consisting of characters '0' and '1' such L[i] == '1' iff (k+1)-mer R[i..i+k] is found in the index.

### Construction

```
Usage:
  bufboss_build [OPTION...]

  -k arg               Order of the de Bruijn graph. Node are k-mers, edges
                       (k+1)-mers. If building from KMC, don't give this.
                       (default: 0)
  -o, --out arg        Output directory. (default: "")
  -a, --add arg        Path to a fasta-file. Adds all (k+1)-mers of the
                       fasta-file to the index. If building from KMC, don't give
                       this. (default: "")
      --add-files arg  Path to a list of fasta-files, one per line. Adds all
                       (k+1)-mers in all the files to the index. If building
                       from KMC, don't give this. (default: "")
  -d, --KMC arg        Build from KMC database (path to a KMC database). The
                       KMC database consists of two files: xxx.kmc_pre and
                       xxx.kmc_suf. You should give only the xxx part here. The
                       database should be built from canonical k-mers (the
                       default behaviour of KMC) (default: "")
  -r, --revcomp        Include reverse complemented k-mers. If building from
                       KMC, don't give this.
  -c, --rrr            Use rrr compression on bit vectors.
  -h, --help           Print instructions.
  -t, --tempdir arg    Directory for temporary working space. (default: "")
```

### Updating

```
Usage:
  bufboss_update [OPTION...]

  -i, --index arg            The directory of the BOSS index. If not given, a
                             new BOSS is built. (default: "")
  -k arg                     If an input index is not given, a new BOSS is
                             built with this k. Otherwise, this k is ignored.
                             (default: 0)
  -o, --out arg              Output directory. (default: "")
  -a, --add arg              Path to a fasta-file. Adds all (k+1)-mers of the
                             fasta-file to the index. (default: "")
      --add-files arg        Path to a list of fasta-files, one per line.
                             Adds all (k+1)-mers in all the files to the index
                             (default: "")
      --add-before-del       If both additions and deletions are given, the
                             deletions are executed first by default. If you
                             want to execute additions first, give this flag.
  -d, --del arg              Path to a fasta-file. Deletes all (k+1)-mers of
                             the fasta-file from the index. (default: "")
      --del-files arg        Path to a list of fasta-files, one per line.
                             Deletes all (k+1)-mers in all the files from the
                             index (default: "")
  -r, --revcomp              Include reverse complemented k-mers.
  -c, --rrr                  Use rrr compression on bit vectors.
      --end-flush            Flush the buffer at the end before writing to
                             disk.
      --count-dummies        Count the number of dummy nodes after the update
  -b, --buffer-fraction arg  If this fraction is x and boss has n nodes, then
                             the buffer is flushed when it has max(n*x,10000)
                             k-mers. (default: 0.01)
  -h, --help                 Print instructions.

```

### Edge existence queries

```
For every input read R, prints to stdout one line L consisting of characters '0' and '1' such L[i] == '1' iff (k+1)-mer R[i..i+k] is found in the index.
Usage:
  bufboss_query [OPTION...]

  -i, --index arg  Path to the directory of the index. (default: "")
  -o, --out arg    Output file. If not given, prints to stdout. (default: "")
  -r, --revcomp    Search reverse-complemented k-mers also.
  -c, --rrr        This option *must* be given if the index was built with
                   rrr compression.
  -h, --help       Print instructions.
  -q, --query arg  Query FASTA-file (default: "")
```

# Limitations

Currently we support only k less or equal to 31. The input files must be in (multi)fasta format.
