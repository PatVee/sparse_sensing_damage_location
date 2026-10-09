# Executed audit and corrections

## Scope

Recovered the uploaded MATLAB/Python source package. Re-ran Python from the two
original raw ZIP archives and executed the corrected shared MATLAB numerical core
in GNU Octave 8.4.0 against all 64 events. MATLAB itself was unavailable.

## Corrections

1. **Missing-value endpoint indexing:** changed `v(badIdx<good(1))` and the matching
   upper-bound expression to index the actual entries in `badIdx`. The old logical
   mask addressed positions in `v` rather than the intended missing sample indices.
2. **Angle interpolation endpoints:** clamp interpolation queries to sampled angle
   endpoints, matching Python's `numpy.interp`, instead of producing NaN outside
   the sampled range.
3. **Angle bin rounding:** explicitly use round-to-even at exact half-centidegree
   ties, matching NumPy. The former MATLAB rounding rule produced substantial
   cross-language alignment differences on these records.
4. **Single-event bootstrap windows:** explicitly retain trials-by-window shape.
   Row-vector indexing could otherwise cause every trial to reuse one sampled event.
5. **Array dimensions:** explicitly reshape event slices as samples-by-channels;
   validate nonempty gauge sets, score dimensions, finite scores and event indices
   before accessing the crack/post-crack records.
6. **Feature coordinate ordering:** match NumPy's upper-triangle Gram ordering for
   arbitrary channel counts. Current one/two-channel zone scores are unaffected.
7. **Execution and presentation:** add the requested `vanersborg_full_pipeline`
   entry point, use `caxis` for the stated MATLAB release compatibility, and derive
   conceptual figure counts/reference budgets from computed results.
8. **Shared implementation:** the MATLAB workflow calls `vanersborg_core.m`, the
   same routines actually exercised by the Octave audit. Vectorised threshold
   crossing and sum/count bin averaging also reduce Octave execution time.

The earlier error trace referred to a different line numbering/function version.
The exact original failure cannot be attributed uniquely to one defect without
that running file. The corrected reduced-array computation has now been executed
successfully on all 64 events; the seven damage indices are retained, not clipped.

## Verified outcomes

- 64/64 raw numeric records aligned and finite.
- Maximum MATLAB-core/Python aligned-strain difference: 3.97e-10.
- Feature parity checked with 1, 2, 5 and 12 channels.
- Primary QR reference indices and POD cumulative energy match Python.
- Maximum virtual-strain difference: 3.97e-10.
- Maximum score difference scaled by max(1, absolute Python score): below 7e-12.
- Primary detection/correct/joint counts: physical 7/7/7; reduced 2/7/2;
  virtual-assisted 7/7/7.
- All 60 configuration/partition outcome rows match exactly across 20 partitions.
- All 20 partition-specific QR reference sets match exactly.
- Empty-array and invalid-event-index regression checks pass.
- Bootstrap draw dimensions checked for windows 1, 3 and 5.

See the execution logs in `audit/` for the actual run output.

## Limits

The MATLAB top-level source was parsed by Octave, but its MATLAB-specific table,
CSV-import and graphics/export APIs were not executed in native MATLAB. Raw CSV
columns were read by pandas and passed unchanged to the MATLAB alignment core for
the Octave comparison. Python's complete extraction-to-figures path was executed.
The audit checks implementation consistency, not statistical independence or the
broader validity of the experimental design. Bootstrap RNG sequences differ
between languages; exact bootstrap Monte Carlo counts are not asserted equal.
