# Interview questions (do this before Day 3)

The plumbing is built and doesn't depend on these answers. Everything smart that sits on
top — which ACMG codes to prioritize, what "trust" looks like, which ancestry cases matter —
does. Get the workflow from the users directly before writing the agent.

## The one question that does most of the work

> "Walk me through the last VUS you classified. What did you open, in what order, and where
> did you lose time?"

That tells you what to automate and what to skip.

## The three scoping questions

1. **Which variant types dominate the queue** — missense SNVs, splicing, CNVs? Pick the
   biggest single bucket and ignore the rest for now. (MVP currently assumes missense +
   nonsense SNVs; confirm or redirect.)
2. **Where does ancestry actually bite?** Matt cares about diagnostic sensitivity in
   non-European patients. Get one real case where thin gnomAD representation made a call
   hard. That case becomes a test fixture for the Day 5 ancestry-confidence layer.
3. **What would they need to see to trust an automated draft?** That answer is the grounding
   spec — it decides how the sources panel and abstention are presented.

## Matt vs. Bridget — decide who to build for

Ask both. Build for whoever's pain is more acute.

- **Matt (clinical):** classify one clinical variant → ACMG 5-tier call. This is what the
  MVP is currently pointed at.
- **Bridget (research):** if her need is sharper, the same engine tilts to research — rank
  candidate variants in an ASD exome by neurodevelopmental gene relevance and de novo
  status. Same retrieval, same grounding, same eval logic; different output.

## What to record here after the interviews

- The concrete last-VUS walkthrough (ordered list of tools + time sinks).
- The chosen variant bucket.
- The real ancestry-hard case (gene + variant + patient ancestry), to add as a fixture.
- The trust checklist (what must be on screen for them to edit rather than redo).
- Decision: Matt-clinical or Bridget-research, with one sentence of why.
