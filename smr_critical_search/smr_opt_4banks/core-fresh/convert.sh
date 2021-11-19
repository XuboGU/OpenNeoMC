for img in *.ppm
do
    c=${img%.*}
    convert "$c.ppm" "$c.png"
done
