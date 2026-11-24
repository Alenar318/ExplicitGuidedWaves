# Colours of the models in the example figures (Okabe-Ito palette, readable with colour
# vision deficiencies), as hex strings that Plots accepts, so the module does not depend
# on Plots. The same model has the same colour in every figure.

"Colour of each model in the figures of `notebooks/examples/`."
const MODEL_COLORS = (closed_form = "#000000",       # black: the exact reference
                      explicit_fem = "#009E73",      # bluish green: this code (as in the paper)
                      hybrid_bem = "#0072B2",        # blue
                      full_3d_bem = "#CC79A7",       # reddish purple
                      implicit_comsol = "#D55E00")   # vermilion
