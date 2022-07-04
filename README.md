# 2D-LeafFlower-Model
Code for the simulation of the medio-lateral section of leaf and floral primordia models in MorphoDynamX.

After installing MorphoDynamX, the models can be compiled by typing:

make -B

If compiling is successfull, flower model can be run with:

MorphoDynamX '--model' 'Model/Folder/00 General' '--addlibrary' 'usrLibFlower.so' 'Flower.mdxv' --run

And leaf model with:

MorphoDynamX '--model' 'Model/Folder/00 General' '--addlibrary' 'usrLibLeaf.so' 'Leaf.mdxv' --run
