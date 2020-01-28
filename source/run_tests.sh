python3 gen_bifurcation.py

mkdir ./density_videos ./density_graphs ./velocity_videos ./velocity_graphs

var="1"
for f in ./out/*.png; do
    python3 blood_flow.py "$f" video --out ./density_videos/density"$var".mp4 --method density -l 100 -g ./density_graphs/density"$var".png
    python3 blood_flow.py "$f" video --out ./velocity_videos/velocity"$var".mp4 --method velocity -l 100 -g ./velocity_graphs/velocity"$var".png
    var=$((var + 1))
done

rm -r ./out
echo "Tests completed!"
