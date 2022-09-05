import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.NumberFormat;
import javax.imageio.ImageIO;

/**
 * Based off of https://stackoverflow.com/questions/2658663
 */
class ImageManipulation {

    private static int percent;

    // TODO: Spec
    private static void manipulateImage(String inputImageFileName, String outputImageFilename) {
        BufferedImage image;
        try {
            image = ImageIO.read(new File(inputImageFileName));
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }
        image = process(image);
        // TODO: Do something better for this case
        if (image == null) {
            System.out.println("image was null, nothing written.");
            return;
        }
        try {
            ImageIO.write(image, "png", new File(outputImageFilename));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    // TODO: Spec
    private static BufferedImage process(BufferedImage old) {
        // Placeholder values for now
        long population = -1;
        int radius = -1;
        boolean populationAndRadiusInitialized = false;
        final String foundPercentCirclesFilename = 
            "/home/alexmijo/PopulationCircles/foundPercentageCircles.txt";
        // TODO: Modularize into a function
        try(BufferedReader br = new BufferedReader(new FileReader(foundPercentCirclesFilename))) {
            String line = br.readLine();
            while (line != null) {
                String[] foundPercentCircle = line.split(" ");
                if (Integer.parseInt(foundPercentCircle[0]) == percent) {
                    population = (long)Double.parseDouble(foundPercentCircle[4]);
                    radius = Integer.parseInt(foundPercentCircle[1]);
                    populationAndRadiusInitialized = true;
                    break;
                }
                line = br.readLine();
            }
        // TODO: Something better when an exception happens, probably
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            return null;
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
        if (!populationAndRadiusInitialized) {
            System.out.println(
                "Desired percentage circle wasn't in " + foundPercentCirclesFilename);
            return null;
        }
        int bigFontHeight = 120;
        int smallFontHeight = 60;
        int verySmallFontHeight = 35;
        double lineSpacing = 1.3;
        int width = old.getWidth();
        int height = old.getHeight() + (int)(bigFontHeight * lineSpacing);
        BufferedImage img = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g2d = img.createGraphics();
        g2d.setColor(Color.BLACK);
        g2d.fillRect(0, 0, width, height);
        g2d.setPaint(Color.RED);
        g2d.setFont(new Font("Serif", Font.BOLD, bigFontHeight));
        NumberFormat commasFormat = NumberFormat.getInstance();
        commasFormat.setGroupingUsed(true); // this will also round numbers, 3 decimal places
        String title = "The smallest possible circle containing " + percent
                       + "% of the world's population";
        String populationString = commasFormat.format(population) + " people";
        String radiusString = "Radius: " + commasFormat.format(radius) +" km";
        g2d.drawImage(old, 0, height - old.getHeight(), null); // ImageObserver not needed
        g2d.drawString(title, 5, bigFontHeight);
        g2d.setFont(new Font("Serif", Font.BOLD, smallFontHeight));
        g2d.drawString(populationString, 5, smallFontHeight + (int)(bigFontHeight * lineSpacing));
        g2d.drawString(radiusString, 5, smallFontHeight + (int)(smallFontHeight * lineSpacing)
                       + (int)(bigFontHeight * lineSpacing));
        
        g2d.setFont(new Font("Serif", Font.BOLD, verySmallFontHeight));
        g2d.setPaint(Color.WHITE);
        FontMetrics fm = g2d.getFontMetrics();
        String dataSourceString = "https://ghsl.jrc.ec.europa.eu/ghs_pop2019.php";
        int dataSourceX = width - fm.stringWidth(dataSourceString) - 5;
        int dataSourceY = height - (int)(verySmallFontHeight * (lineSpacing - 1)) - 5;
        g2d.drawString(dataSourceString, dataSourceX, dataSourceY);
        String dataResolutionString = "30 arcsecond resolution";
        int dataResolutionX = width - fm.stringWidth(dataResolutionString) - 5;
        int dataResolutionY = dataSourceY - (int)(verySmallFontHeight * lineSpacing);
        g2d.drawString(dataResolutionString, dataResolutionX, dataResolutionY);
        String dataYearString = "2015 population data";
        int dataYearX = width - fm.stringWidth(dataYearString) - 5;
        int dataYearY = dataResolutionY - (int)(verySmallFontHeight * lineSpacing);
        g2d.drawString(dataYearString, dataYearX, dataYearY);
        g2d.dispose();
        return img;
    }

    public static void main(String[] args) {
        for (percent = 90; percent <= 90; percent++) {
            String inputImageFileName =
                "/mnt/c/Users/Administrator/Desktop/linuxPythonMadePercentMaps/" + percent + 
                "PercentCircle.png";
            String outputImageFileName =
                "/mnt/c/Users/Administrator/Desktop/linuxPythonMadePercentMaps/" + percent +
                "PercentCircleWithText.png";
            manipulateImage(inputImageFileName, outputImageFileName);
        }
    }
}