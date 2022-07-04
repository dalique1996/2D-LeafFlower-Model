# 2D-LeafFlower-Model
Code for the simulation of the medio-lateral section of leaf and floral primordia models in MorphoDynamX.

After installing MorphoDynamX, the models can be compiled by typing:

make -B

If compiling is successfull, the models can be run with:

MorphoDynamX '--model' 'Model/Root/01 Root' '--addlibrary' 'usrLibRoot.so' Root_model.mdxv
