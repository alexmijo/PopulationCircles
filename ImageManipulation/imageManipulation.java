import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;

import javax.imageio.ImageIO;
import javax.swing.JPanel;

/**
 * @see https://stackoverflow.com/questions/2658663
 */
class ImageManipulation extends JPanel {

    private static final int percent = 2;
    private static final int population = 147265836;
    private static final int radius = 177;

    private static void manipulateImage(String inputImageFileName, String outputImageFilename) {
        BufferedImage image;
        try {
            image = ImageIO.read(new File(inputImageFileName));
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }
        image = process(image);
        try {
            ImageIO.write(image, "png", new File(outputImageFilename));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static BufferedImage process(BufferedImage old) {
        int bigFontHeight = 120;
        int smallFontHeight = 60;
        double lineSpacing = 1.3;
        int w = old.getWidth();
        int h = old.getHeight() + (int)(bigFontHeight * lineSpacing);
        BufferedImage img = new BufferedImage(w, h, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g2d = img.createGraphics();
        g2d.setColor(Color.BLACK);
        g2d.fillRect(0, 0, w, h);
        g2d.setPaint(Color.RED);
        g2d.setFont(new Font("Serif", Font.BOLD, bigFontHeight));
        NumberFormat commasFormat = NumberFormat.getInstance();
        commasFormat.setGroupingUsed(true); // this will also round numbers, 3 decimal places
        String title = "The smallest possible circle containing " + percent
                       + "% of the world's population";
        String populationString = commasFormat.format(population) + " people";
        String radiusString = "Radius: " + commasFormat.format(radius) +" km";
        g2d.drawImage(old, 0, h - old.getHeight(), null); // ImageObserver not needed
        g2d.drawString(title, 5, bigFontHeight);
        g2d.setFont(new Font("Serif", Font.BOLD, smallFontHeight));
        g2d.drawString(populationString, 5, smallFontHeight + (int)(bigFontHeight * lineSpacing));
        g2d.drawString(radiusString, 5, smallFontHeight + (int)(smallFontHeight * lineSpacing)
                       + (int)(bigFontHeight * lineSpacing));
        g2d.dispose();
        return img;
    }

    public static void main(String[] args) {
        String inputImageFileName = "/mnt/c/Users/Administrator/Desktop/linuxPythonMadePercentMaps/"
                                    + percent + "PercentCircle.png";
        String outputImageFileName =
            "/mnt/c/Users/Administrator/Desktop/linuxPythonMadePercentMaps/" + percent
            + "PercentCircleWithText.png";
        manipulateImage(inputImageFileName, outputImageFileName);
    }
}