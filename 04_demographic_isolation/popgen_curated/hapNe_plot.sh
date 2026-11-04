#!/usr/bin/env bash


# Summarize population variation and demographic history.
# Set input, reference and tool paths for your environment before execution.

cd /path/to/hapNe || exit 1

python3 hapNeIBD.plot.multi.py \
-i 4panel_HapNeIBD/25Malay/HapNe/hapne.csv \
4panel_HapNeIBD/36Bajau/HapNe/hapne.csv \
4panel_HapNeIBD/EGDP.Borneo/HapNe/hapne.csv \
4panel_HapNeIBD/EGDP.Luzon/HapNe/hapne.csv \
4panel_HapNeIBD/EGDP.Myanmar/HapNe/hapne.csv \
4panel_HapNeIBD/EGDP.Vietnam/HapNe/hapne.csv \
4panel_HapNeIBD/Pacific.Barito/HapNe/hapne.csv \
4panel_HapNeIBD/Pacific.Besemah/HapNe/hapne.csv \
4panel_HapNeIBD/Pacific.Kalimantan/HapNe/hapne.csv \
4panel_HapNeIBD/Pacific.Mamanwa/HapNe/hapne.csv \
4panel_HapNeIBD/Pacific.Manobo/HapNe/hapne.csv \
4panel_HapNeIBD/Pacific.Semende/HapNe/hapne.csv \
4panel_HapNeIBD/PanMLS.Jakun/HapNe/hapne.csv \
4panel_HapNeIBD/PanMLS.Temiar/HapNe/hapne.csv \
 -o 4panel_HapNeIBD.plot.pdf --min_time 1 --max_time 100 \
 --figsize 10 7

python3 hapNeIBD.plot.multi.py \
-i 4panel_HapNeIBD/36Bajau/HapNe/hapne.csv \
4panel_HapNeIBD/Pacific.Barito/HapNe/hapne.csv \
 -o 4panel_HapNeIBD.Pacific.Barito_36Bajau.plot.pdf --min_time 1 --max_time 100 \
 --figsize 10 7
