# Speaker script — Structural Clustering of Cancer-Risk Gene Variants

*Matches the 12-slide deck. Read it over a few times rather than word-for-word.*

---

### Slide 1 — Title

Good [morning/afternoon], thanks for having me. Today I'll walk through a
project I've been building: using protein structure to help prioritize
genetic variants that doctors can't currently classify — combining
AlphaFold's predicted structures with a bit of machine learning. I'll cover
what we built, why, and what we found.

---

### Slide 2 — The unclassified-variant problem

Genetic testing often turns up variants that aren't clearly disease-causing
or harmless — these are Variants of Uncertain Significance, or VUS.
Classifying them formally requires new clinical or lab evidence;
computational evidence like this is the weakest tier and can't confirm
anything alone. So to be upfront: this project generates leads worth
investigating — it doesn't diagnose anyone.

---

### Slide 3 — A two-phase approach

The project has two phases. Phase one is geometric: flag variants close to
a known pathogenic spot in 3D space, even if they're far apart in the
linear sequence. Phase two asks whether that signal holds up
statistically — pooling it across many genes and testing whether it
predicts pathogenicity on genes it's never seen.

---

### Slide 4 — Data sources and scale

For each gene we needed three things: ClinVar for known classifications,
AlphaFold for 3D structure, UniProt for sequence and domain info. We ended
up with 24 genes, screened from 125 candidates, covering about 48,000
classified variants — after checking each gene for numbering errors and
cross-gene contamination, which a few raw datasets actually had.

---

### Slide 5 — Flagging by structure, not sequence

Core rule for phase one: flag a VUS if it's within 6 angstroms of a known
pathogenic residue in 3D, but more than 10 positions away in the linear
sequence. That sequence-distance condition matters — without it, almost
everything flagged was just next to a pathogenic spot in the sequence
anyway, which isn't interesting. Requiring both is what isolates cases the
sequence alone would never suggest.

---

### Slide 6 — What the model actually checks

Phase two reduces every variant to seven numbers: 3D distance and sequence
distance to a pathogenic residue, their ratio and difference, how many
pathogenic residues are nearby, structural confidence, and whether the
region is disordered. We train on variants we already know are pathogenic
or benign, then score new ones. We used simple models — logistic
regression and a small random forest — specifically so we could see which
measurements actually mattered, instead of trusting a black box.

---

### Slide 7 — Cross-gene generalization

We tested generalization by holding out one gene at a time, training on
the rest, and checking performance on the held-out gene — repeated for all
24. Average balanced accuracy: about 0.72, where 0.50 is chance. But it's
uneven — genes like SCN2A, TP53, and BRCA1 hit 0.8 or higher; others like
RB1, COL5A1, and MYBPC3 land near chance. That spread is the honest
finding — this doesn't work the same way on every gene.

---

### Slide 8 — Ruling out a confound

Before trusting this, we had to rule something out: maybe the model wasn't
detecting real 3D clustering — just learning "this is an important,
well-folded region," a weaker signal. So we removed the
structural-confidence measurement entirely and reran everything. Accuracy
barely moved — 0.72 to 0.71. That tells us the clustering signal is doing
real, independent work — not just riding on a proxy for structure quality.

---

### Slide 9 — Benchmark against AlphaMissense

I benchmarked against AlphaMissense, DeepMind's genome-scale predictor.
Here's why: my own accuracy number means little in isolation — I needed an
independent reference to know if it's meaningful. AlphaMissense is a good
one because it's built from completely different evidence and knows
nothing about my 3D-clustering approach — so if it agrees with my
candidates more than chance, that's real outside validation, not just my
model overfitting itself.

We're not trying to do the same job, though: AlphaMissense is broad and
general, trained on tens of millions of variants; mine only checks one
thing — 3D closeness to a known pathogenic spot. So it's no surprise it's
more accurate alone, 0.86 versus my 0.72.

Correlation between us was moderate — expected, since we use different
evidence and should agree on obvious cases, diverge on subtle ones. The
number I care about most: of my highest-confidence candidates,
AlphaMissense independently agrees on about half — versus a baseline of
about a third across all uncertain variants. So my flagging is
meaningfully enriched for cases an independent model also flags as
concerning.

The disagreements are interesting too — either new signal, or false
positives worth checking by hand. I don't know which yet, and I'd rather
say that plainly than oversell it. But the agreement rate is real evidence
this structural signal adds something, not just echoes a bigger model.

---

### Slide 10 — Limitations

Some honest limitations: this dataset is small next to genome-scale
tools — tens of thousands of variants versus tens of millions. Performance
varies a lot gene to gene, so no single number fairly represents it. We're
using static structures, so we miss protein flexibility and complex
behavior. And again: this is a hypothesis-generation tool — every output
is a candidate, never a classification.

---

### Slide 11 — Where this goes next

Right now I'm manually following up on some flagged candidates — checking
population frequency, other predictors, and existing literature. Next
steps: more genes for a firmer statistical picture, and richer structural
features beyond distance and confidence. But the core point stays the
same — this is a lead-generation tool. Its value is in surfacing
candidates for someone to actually investigate.

---

### Slide 12 — Thank you, Dr. Pellegrini

Thank you, Dr. Pellegrini, for hearing this out — and for the mentorship
that got me into this kind of work in the first place. Happy to take
questions.
