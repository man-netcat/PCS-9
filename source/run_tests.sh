var="1"
for f in ./out/*.png; do
    python3 blood_flow_video_test.py "$f" "$var.mp4"
    var=$((var + 1))
done
