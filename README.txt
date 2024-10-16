1.	MT_Toruqe_Jump_Pol_II_Position.m
Script to load the converted data including time, magnet angle and extension from magnetic tweezers instrument for torque jump experiments described in figure 2, converting the extension signal into Pol II position on the template.

2.	MT_Torque_Jump_Pause_free_velocity.m
Script to load the converted data including time, magnet angle and extension from magnetic tweezers instrument for torque jump experiments described in figure 2, analyzing trace by trace and getting the pause-free-velocity of Pol II.


3.	MT_Torque_Jump_Active_Fraction.m
Script to load the converted data including time, magnet angle, magnet height and extension from magnetic tweezers instrument for torque jump experiments described in figure 2, analyzing trace by trace and getting the active fraction of Pol II.

4.	MT_Pol_II_Position_Through_Single_Nucleosome.m
Script to load the converted data including time, magnet angle and extension from magnetic tweezers instrument for PoI II transcription on single nucleosomal template experiments described in figure 3, analyzing trace by trace and converting the extension signal into Pol II position on the single nucleosomal template.

5.	MT_Pol_II_Position_Through_Nucleosome_Array.m
Script to load the converted data including time, magnet angle and extension from magnetic tweezers instrument for PoI II transcription on nucleosome array template experiments described in figure 3, analyzing trace by trace and converting the extension signal into Pol II position on the nucleosome array template.

6.	MT_Pol_II_Dwell_Time_Histogram.m
Script to load the converted Pol II position v.s. time data (on single nucleosome or nucleosome array template), analyzing the dwell time pattern of Pol II inside the NPE as in figure 3.

7.	MT_Average_Trace_Plotting.m
Script to load the converted Pol II position v.s. time data (on single nucleosome or nucleosome array template), analyzing the averaged position of Pol II on the template as in figures 4 and 5.

8.	MT_Pol_II_Topo_Through_Nucleosome.m
Script to load the converted data including time, magnet angle and extension from magnetic tweezers instrument for PoI II transcription on nucleosome array with topo II experiments described in figure 5, analyzing trace by trace and comparing the initial and the final extension-turns curves, getting the number of nucleosome Pol II transcribed through.

9.	Nucleosome_array_hat_curve.m
A function, with the input as the coordinates of 4 specific positions of the extension-turns curve of a nucleosome array, and returns the shape of such extension-turns curve. Used in “MT_Pol_II_Topo_Through_Nucleosome.m”.

10.	Hat_curve_after_transcribing_through_ita_nucleosomes.m
A function, returns the shape of a predicted extension-turns curve based on an initial extension-turns curve after Pol II transcribed through a certain number (ita) of nucleosomes. Used in “MT_Pol_II_Topo_Through_Nucleosome.m”.


11.	MT_Pol_II_Position_Through_Nucleosome_Array_Continuous_Winding.m
Script to load the converted data including time, magnet angle and extension from magnetic tweezers instrument for PoI II transcription on nucleosome array template with continuous winding described in figure 5, analyzing trace by trace and converting the extension signal into Pol II position on the nucleosome array template.

12.	AOT_transcription_data_analysis.m
Script to process the AOT data and determine the torque and position during transcription by Pol II.  Saves the maximum torque Pol II can generate and the position versus time traces.

13.	pol_II_position_vs_time_analysis.m
Script to identify pauses in the pol II versus time traces.

14.	calculate_velocity_position.m
Function that calculates the smoothed velocity and position versus time from the raw position versus time trace.

15.	AOT_animate_pol_II_data.m
Script to make an animation of the pol II rotations versus time data.
