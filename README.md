# H1 3DBMPP refactor

This is a refactor of the  3D Beta-barrel Membrane Protein Predictor code  by Tian, Wei and Lin, Meishan and Tang, Ke and Liang, Jie and Naveed, Hammad
The source code is available at : http://sts.bioe.uic.edu/3dbmpp/


# H3 Full Citation
<cite>
@article {Tian201716817,
	author = {Tian, Wei and Lin, Meishan and Tang, Ke and Liang, Jie and Naveed, Hammad},
	title = {High-resolution structure prediction of β-barrel membrane proteins},
	year = {2018},
	doi = {10.1073/pnas.1716817115},
	publisher = {National Academy of Sciences},
	abstract = {β-Barrel membrane proteins (βMPs) are drawing increasing attention because of their promising potential in bionanotechnology. However, their structures are notoriously hard to determine experimentally. Here we develop a method to achieve accurate prediction of βMP structures, including those for which no prediction has been attempted before. The method is general and can be applied to genome-wide structural prediction of βMPs, which will enable research into bionanotechnology and drugability of βMPs.β-Barrel membrane proteins (βMPs) play important roles, but knowledge of their structures is limited. We have developed a method to predict their 3D structures. We predict strand registers and construct transmembrane (TM) domains of βMPs accurately, including proteins for which no prediction has been attempted before. Our method also accurately predicts structures from protein families with a limited number of sequences and proteins with novel folds. An average main-chain rmsd of 3.48 {\r A} is achieved between predicted and experimentally resolved structures of TM domains, which is a significant improvement (\&gt;3 {\r A}) over a recent study. For βMPs with NMR structures, the deviation between predictions and experimentally solved structures is similar to the difference among the NMR structures, indicating excellent prediction accuracy. Moreover, we can now accurately model the extended β-barrels and loops in non-TM domains, increasing the overall coverage of structure prediction by \&gt;30\%. Our method is general and can be applied to genome-wide structural prediction of βMPs.},
	issn = {0027-8424},
	URL = {http://www.pnas.org/content/early/2018/01/25/1716817115},
	eprint = {http://www.pnas.org/content/early/2018/01/25/1716817115.full.pdf},
	journal = {Proceedings of the National Academy of Sciences}
}
</cite>
#H2 Requirements
Sidechain prediction relies on [SCWRL](http://dunbrack.fccc.edu/scwrl4/). To download scwr, you must apply for a license on their website.

To install python requirements,

#H2 Usage

To predict beta-barrel structure, you must provide a fasta file with the protein sequence as well as the identified transmembrane strands. (See the example_fasta.txt file for reference).

We recommend using [PRED-TMBB2](https://bio.tools/PRED-TMBB2) to predict the transmembrane strands.

Than call:
`python main.py -f [fastfile]`
to predict the beta barrel structure.  the new pdb files can be found in `/output`


