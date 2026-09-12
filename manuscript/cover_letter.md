Dulian Zeqiraj
Department of Energy Resources
Faculty of Geology and Mining
Polytechnic University of Tirana
Rruga e Elbasanit, Tirana 1001, Albania
dulian.zeqiraj@fgjm.edu.al

The Editors
Journal of Hydrology

Dear Editors,

I am submitting for your consideration the manuscript "What is a nitrate measurement worth? Fisher information and goal-oriented monitoring network design in alluvial aquifers".

The paper asks a question that monitoring practice usually answers by assertion: what is the data a network already holds actually worth, and where should the next well go. It answers it for one aquifer, the Fushë-Kuqe alluvial system in northwestern Albania, using the whole of the available field record and releasing all of it.

Three things make it worth your readers' time, and I would rather state them plainly than dress them up.

The first result is negative and it corrects my own earlier work. On the 24 nitrate wells inside this aquifer, no method I could test predicts a withheld well better than the network mean: not inverse-distance weighting, not ordinary kriging at any correlation length between 1 and 16 km, not linear regression on the available covariates, and not a transient transport model with a loading field fitted to the data. Two recent papers of mine on this same dataset report a nitrate prediction skill of R² = 0.868. That is an in-sample number. Cross-validated, nothing on these data comes close, and I think it matters that this is said in the literature rather than left for a reader to discover. The same computation shows why: under a uniform regional flux, effective porosity cancels out of the steady concentration field, so the snapshot carries a very small fraction of the Fisher information about porosity, and the corollary drawn in that earlier work, that one concentration observation is worth about three pumping tests, comes out inverted when it is computed rather than argued.

The second is a small theorem with a practical edge. Several proposals in the vulnerability-assessment literature weight a D-optimal design criterion by a vulnerability field in order to steer the design toward vulnerable ground. That cannot work: D-optimality is invariant to any fixed reweighting of the parameter space, and the proof is two lines. What does work is goal orientation, targeting the variance of the vulnerability-weighted contaminant arrival time, and on these data it reduces that variance by a large multiple of what the classical criterion achieves.

The third is a reproducibility commitment I would ask you to hold me to. Every number in the paper is emitted by one run of one pipeline into a single results file, and the manuscript is assembled by a script that refuses to compile if any number in the text lacks a value in that file. The code, the data in flat form and the results file are released on GitHub and archived on Zenodo, and the MATLAB listings are printed in full after the references. The conclusions are replicated on an independent public USGS dataset from the Mississippi Alluvial Plain.

The manuscript is original, is not under consideration elsewhere, and has a single author with no competing interests. It overlaps two of my published papers only in the sense that it takes their porosity fusion as a prior input and cites them as such; the fusion itself is not re-derived here, and where the present computation disagrees with them I say so and show the calculation.

Yours sincerely,

Dulian Zeqiraj
