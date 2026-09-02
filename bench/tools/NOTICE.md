`emps.c` is the standard MPS decompressor from the Netlib LP archive
(https://www.netlib.org/lp/data/emps.c), fetched verbatim. It only unpacks
the Netlib set's compact on-disk encoding back into plain-text MPS; it
performs no optimization and is not a solver or solver component. Per
`BUILD_PLAN_V2.md`'s "not built upon any existing open source solver
library" rule, this is a data-format tool, not a dependency of the solver
itself, and is never linked into `inferno_core` or `solver`.
