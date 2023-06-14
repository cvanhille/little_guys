README FILE EXPLAINING MODIFICATIONS IMPLEMENTED ON 09/06/2023 TO:
  fix_bond_create.h
  fix_bond_create.cpp

Added the inter_mol flag option to the fix. This flag does not affect any other flags. Its only purpose is to specify that bond creation may only occur between atoms belonging to different molecule: their molecule identifiers need to be different. All other constraints should still hold. This means that a few lines were added to the FixBondCreate::post_integrate() function in fix_bond_create.cpp to skip pairs with identical molecular identifiers.