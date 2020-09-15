# RaxModeler
Fast selection of models of evolution using RAxML for ML/Bayesian phylogenies. This is a reverse-engineered version of [ProtTest2](https://academic.oup.com/bioinformatics/article/21/9/2104/408994) (Java). Rather than run the single-threaded much slower PhyML (Java) for model testing, this Python wrapper ran multiple instances of RaXML in the background (C++).   
