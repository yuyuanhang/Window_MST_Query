# Window MST
## Introduction
This project includes:

(1) four query algorithms for solving window-MST query in temporal graphs.

- an online algorithm ```Online```.
- a $\mathrm{\Gamma}$ index-based query algorithm ```Query```.
- a $\mathrm{\Gamma}_d$ index-based query algorithm ```QueryD```.
- a $\mathrm{\Gamma}_b$ index-based query algorithm ```QueryB```.

(2) four index construction algorithms.

- an index construction algorithm ```ERTSConstruction``` for the index $\mathrm{\Gamma}$.
- an optimized index construction algorithm ```ERTSConstruction*``` for the index $\mathrm{\Gamma}$.
- an index construction algorithm ```DIndexConstruction``` for the index $\mathrm{\Gamma}_d$.
- an index construction algorithm ```BIndexConstruction``` for the index $\mathrm{\Gamma}_b$.

(3) an multi-threaded index construction algorithm for the index $\mathrm{\Gamma}$.

(4) an index maintenance algorithm for the index $\mathrm{\Gamma}_b$ under incremental updates.

(5) a window-MST query generator.

## Usage
### Compile program
```zsh
mkdir build
cd build
cmake ..
make
```

### Query Generation
To invoke the window-MST query generator, the following command line is required:
```zsh
./build/MST -gq [data graph] [window_ratio] [num_queries] [query]
```
- The parameter ```[data graph]``` is the path of dataset. Specifically, the dataset is stored in a binary file: (1) the first 4 bytes (an int) indicate the number of edges. (2) an array of Edge structs follows, written in binary.
- The parameter ```[window_ratio]``` is the ratio of query time window size to $t_{max}$, which is a floating-point number in ```(0, 1]```.
- The parameter ```[num_queries]``` specifies the number of queries.
- The parameter ```[query]``` specifies the prefix of the storage path of generated queries.

For example, the following command line is used to generate 1,000 window-MST queries with size $0.2t_{max}$ on the ```SocBitcoin``` dataset:

```zsh
./build/MST -gp datasets/socbitcoin.bin 0.2 1000 query/socbitcoin
```

### Construction
To invoke the index construction algorithms, the following command line is required:
```zsh
./build/MST -bi [data graph] [method_type] [index]
```
- The parameter ```[data graph]``` is the path of dataset. Specifically, the dataset is stored in a binary file: (1) the first 4 bytes (an int) indicate the number of edges. (2) an array of Edge structs follows, written in binary.

- The parameter ```[method_type]``` specifies the type of index. Specifically, (1) ```base``` for ```ERTSConstruction```, (2) ```op``` for ```ERTSConstruction*```, (3) ```d``` for ```DIndexConstruction```, and (4) ```b``` for ```BIndexConstruction```.

- The parameter ```[ndex]``` specifies the prefix of the storage path of index.

For example, the following command line is used to build the index $\mathrm{\Gamma}_b$ on the ```SocBitcoin``` dataset:

```zsh
./build/MST -bi datasets/socbitcoin.bin b idx/socbitcoin
```

To invoke the multi-threaded index construction for the index $\mathrm{\Gamma}$, the following command line is required:
```zsh
./build/MST -bip [data graph] [num_thread] [index]
```
Here, parameters are the same as described above, the only difference is that the parameter ```[num_thread]``` is the number of threads.

To invoke the index maintenance algorithm for the index $\mathrm{\Gamma}_b$ under incremental updates, the following command line is required:
```zsh
./build/MST -bii [data graph] [edge_ratio]
```
Here, parameters are the same as described above, the only difference is that the parameter ```[edge_ratio]``` is the initial index span, which is a floating-point number in ```(0, 1]```.

### Search
To invoke query algorithms, the following command line is required:
```zsh
./build/MST -eq [data graph] [query] [method_type] [index]
```
- Similarly, the parameter ```[data graph]``` is the path of dataset.

- The parameter ```[query]``` is the path of query.

- The parameter ```[method_type]``` specifies the type of algorithm. Specifically, (1) ```online``` for ```Online```, (2) ```base``` for ```Query```, (3) ```d``` for ```QueryD```, and (4) ```b``` for ```QueryB```.

- The parameter ```[index]``` is the prefix of the storage path of index. Note that if the online algorithm is used, this parameter can be any string.

For example, the following command line is used to perform ```QueryB``` on the ```SocBitcoin``` dataset to resolve the queries generated above:

```zsh
./build/MST -eq datasets/socbitcoin.bin query/socbitcoin_0.2.bin b idx/socbitcoin
```

## Datasets
The information of used real-world datasets is provided in our paper.
