
############################################################################
# "Steady" examples

create_macro_ebacoli_executable(driver95-simple steady 1 0 0)
create_macro_ebacoli_executable(driver95-curve steady 1 0 0)
create_macro_ebacoli_executable(driver95-trimesh steady 1 0 0)

############################################################################
# Cahn--Allen equation examples

create_macro_ebacoli_executable(driver95-simple cahn_allen 1 0 0)
create_macro_ebacoli_executable(driver95-curve cahn_allen 1 0 0)
create_macro_ebacoli_executable(driver95-trimesh cahn_allen 1 0 0)

############################################################################
# Reaction-diffusion system examples

create_macro_ebacoli_executable(driver95-simple rcdsys 4 0 0)
create_macro_ebacoli_executable(driver95-curve rcdsys 4 0 0)
create_macro_ebacoli_executable(driver95-trimesh rcdsys 4 0 0)

############################################################################
