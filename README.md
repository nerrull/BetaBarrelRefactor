# 3D-BMPP refactor

This is a refactor of the  3D Beta-barrel Membrane Protein Predictor code by Will Tian. 
The original source code is available at : http://sts.bioe.uic.edu/3dbmpp/

### Full Citation
```
    @article {Tian201716817,
        author = {Tian, Wei and Lin, Meishan and Tang, Ke and Liang, Jie and Naveed, Hammad},
        title = {High-resolution structure prediction of β-barrel membrane proteins},
        year = {2018},
        doi = {10.1073/pnas.1716817115},
        publisher = {National Academy of Sciences},
        issn = {0027-8424},
        URL = {http://www.pnas.org/content/early/2018/01/25/1716817115},
        eprint = {http://www.pnas.org/content/early/2018/01/25/1716817115.full.pdf},
        journal = {Proceedings of the National Academy of Sciences}
    }
```

## Requirements
Sidechain prediction relies on [SCWRL](http://dunbrack.fccc.edu/scwrl4/). To download scwrl, you must apply for a license on their website. 
We then recommend unpacking it into a folder name `scwrl`in the top level folder of this repo. You can also change the path in `defs.py`.

Python requirements `numpy, pandas, biopython`

To install required python packages: `pip install -r requirements.txt`

## Compiling
Part of the application is written in C and needs to be compiled.

```
cd bb_register/pred_combine/
make
```

If you're getting an error that says something like `Can't find file coonstruction_inputs/resuls/regs/*.regs` you probably forgot to do this step 


## Usage

To predict beta-barrel structure, you must provide a fasta file with the protein sequence as well as the identified transmembrane strands. 
See the `example_fasta.txt` file for reference.
You can also provide a fasta file with multiple sequences and they will be processed individually.

We recommend using [PRED-TMBB2](https://bio.tools/PRED-TMBB2) to predict the transmembrane strands.

### Running the code 

From a shell call: `python main.py -f [fastfile]`

to predict the beta barrel structure.  The newly created pdb files can be found in `output/` and can be inspected using [PyMOL](https://pymol.org/2/). 

#### Other parameters:
`--level` or `-l` : \[1-5\] or -1 
 *This changes the set of hyperparameters used during register prediction.

A detailed description of these parameters can be found on page 2 of the appendix in `docs/`.
If the parameter is set to -1, a predicted structure is generated for each level. The level of parameters used is identified in output files by a `_l0#` suffix.

When the parameter is not defined, we infer the level based on the number of transmembrane strands (N) in the provided input file. 
* For N<16, 2 barrels are predicted using both levels 1 and 2
* For 16 <= N < 20, 2 barrels are predicted using levels 3 and 4

Quick reference :

* 1 : Small BMPs (N < 16) without inplugs or outclamps
* 2 : Small BMPs (N < 16) with inplugs or outclamps
* 3 : Medium oligomeric BMPs (16 <= N < 20)
* 4 : Medium monomeric BMPs (16 <= N < 20)
* 5 : Large BMPs (N >= 20)




