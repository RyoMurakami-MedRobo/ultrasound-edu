# HW3: DAS experiments and execution replay

The simulator supports the DAS learning experiments in Question 4 of
HW3_2025. It uses MUST/SIMUS or analytic mock RF in MATLAB, and analytic mock
RF in Python. It does **not** run Field II or include the linked kidney data,
`DAS.m`, `imageGen_256`, or `sim_kidney_linear_array_256.m` from the assignment.
Those files were not supplied. Results here are simulator experiments, not
reproductions of the original Field II kidney exercise. The DMAS bonus and
curvilinear-array exercise are not implemented by this change.

## Run and see the computation

1. Start `main_gui` in MATLAB, or `python -m ultrasound_das` from `python/`.
2. Choose a point phantom and press **[1] Simulate**. The Wave animation tab
   opens and plays one sweep automatically. Play replays it; Stop and the
   time slider allow inspection. This wave is a propagation illustration,
   not an animation of a beamforming algorithm.
3. Copy the updated custom template, rename the function, and edit its DAS
   loop. In Detailed mode set **Custom function name** (MATLAB) or
   **Custom module:function** (Python). Choose Custom DAS as B.
4. Press **[2] Beamform / Compare**. A and B are reconstructed, then the
   **DAS execution replay** tab opens and plays B automatically. Choose A
   or B and press **Replay DAS** to compare implementations. If B is None,
   A is selected automatically for replay.
5. Select a reconstruction point in B-mode or the delay tab and replay.
   The replay executes the selected implementation over the entire depth
   line at that x coordinate, using the selected transmit event. For a
   multi-focus image it shows that individual event, not the zoned composite.

The **Delay curve & alignment** tab places B-mode A beside the acquired RF.
Click or drag the orange circle on B-mode to choose a reconstruction pixel.
The orange curve on RF shows the per-element samples used for that same
pixel. Before/after alignment use the same orientation as the RF panel:
receive element on x, time on y, with time and depth increasing downward.
After alignment, the selected samples meet at relative time zero and the red
trace at the right is their coherent sum.

The four panels show acquired RF with the implementation's sampling delays,
the actual delayed and weighted channel contributions revealed in order,
the signed accumulating sum along depth, and the final linear envelope
returned by the implementation. Signals of the same sign reinforce and
opposite signs cancel. The final envelope is labeled as a final result; it
is not a partial Hilbert envelope. Replay is a recorded execution of the
scan line, after reconstruction, not a live progress indicator for the
whole image.

## Keep the execution hooks in your implementation

MATLAB accepts the unchanged seven inputs and an optional third output:
`[image, delays, trace]`. Python adds `want_trace=True` and returns the same
three values. Normal calls retain their previous return values and images.

Keep the trace capture immediately around the **actual** accumulator update:

```matlab
previous = bf;
bf = bf + w .* val; % your implementation
if nargout > 2
    trace.contributions(:, e) = bf - previous;
end
```

```python
previous = bf.copy() if want_trace else None
bf += w * val  # your implementation
if want_trace:
    contributions[:, e] = bf - previous
```

`trace.contributions` is `[Npix, Nel]`, column-major pixel order;
`trace.coherent` is the signed sum reshaped to the image dimensions.
Keep the delay output consistent with your sampling, using NaN for excluded
channels. The replay verifies that channel increments sum to the coherent
result. Older custom functions still reconstruct but report **Trace
unavailable** in replay until these hooks are added. Replay never substitutes
a reference beamformer for a missing custom trace.

## Assignment experiment map

| Exercise | Simulator procedure |
| --- | --- |
| DAS geometry | Edit transmit time, receive distance and RF sample selection in the custom template; compare A/B and replay both. |
| 128 / 256 / 512 elements | Select Elements, simulate again, then compare. Pitch stays fixed, so physical aperture changes. Set receive f-number to 0 for a full-aperture comparison; a dynamic aperture can mask the effect of more elements. |
| Sampling frequency | Change fs/fc in Detailed mode and simulate again; compare nearest-neighbor and linear sampling. |
| Hilbert transform | Inspect the signed sum and final envelope panels; edit the envelope step and compare. |
| Point resolution | Use a single point and the Comparison tab's lateral -6 dB FWHM and profile. Increase reconstruction grid density for finer measurements. |
| Apodization | Select rect/hann for Reference DAS. The starter Custom DAS intentionally uses rectangular weights; implement the weighting in your copy to compare it. |

Preserve the same phantom, transmit mode, pitch, frequency and reconstruction
grid across comparisons. Save screenshots of B-mode and Comparison for the
image and resolution experiments. Full Field II assignment compatibility
requires the original exercise files and an explicit RF/time-origin adapter.
