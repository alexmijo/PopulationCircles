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
import java.text.DecimalFormat;
import java.text.NumberFormat;
import javax.imageio.ImageIO;
import static java.util.Map.entry;
import java.util.Map;
import java.lang.Math;


/**
 * Based off of https://stackoverflow.com/questions/2658663
 */
class ImageManipulation {

    private static String percent;
    private static int circleNum;

    private static Map<String, String> countries = Map.ofEntries(
        entry("89.0", ""),
        entry("89.1", ""),
        entry("89.2", ""),
        entry("89.3", ""),
        entry("89.4", ""),
        entry("89.5", ""),
        entry("89.6", ""),
        entry("89.7", ""),
        entry("89.8", ""),
        entry("89.9", "Sweden"),
        entry("90.0", "Sweden"),
        entry("90.1", "Sweden"),
        entry("90.2", "Sweden"),
        entry("90.3", "Sweden"),
        entry("90.4", "Sweden"),
        entry("90.5", "Sweden"),
        entry("90.6", "Sweden"),
        entry("90.7", "Sweden"),
        entry("90.8", "Sweden"),
        entry("90.9", "Sweden"),
        entry("91.0", "Sweden"),
        entry("91.1", "Sweden"),
        entry("91.2", "Germany"),
        entry("91.3", "Germany"),
        entry("91.4", "Germany"),
        entry("91.5", "Germany"),
        entry("91.6", "Germany"),
        entry("91.7", "Germany"),
        entry("91.8", "Germany"),
        entry("91.9", "Germany"),
        entry("92.0", "Germany"),
        entry("92.1", "Germany"),
        entry("92.2", "Germany"),
        entry("92.3", "Switzerland"),
        entry("92.4", "Switzerland"),
        entry("92.5", "France"),
        entry("92.6", "Germany"),
        entry("92.7", "Germany"),
        entry("92.8", "Czechia"),
        entry("92.9", "Czechia"),
        entry("93.0", "Czechia"),
        entry("93.1", "Czechia"),
        entry("93.2", "Germany"),
        entry("93.3", "Germany"),
        entry("93.4", "the Arctic Ocean"),
        entry("93.5", "the Arctic Ocean"),
        entry("93.6", "the Arctic Ocean"),
        entry("93.7", "the Arctic Ocean"),
        entry("93.8", "the Arctic Ocean"),
        entry("93.9", "the Arctic Ocean"),
        entry("94.0", "the Arctic Ocean"),
        entry("94.1", "the Arctic Ocean"),
        entry("94.2", "Russia"),
        entry("94.3", "Russia"),
        entry("94.4", "Russia"),
        entry("94.5", "Russia"),
        entry("94.6", "Russia"),
        entry("94.7", "Russia"),
        entry("94.8", "Russia"),
        entry("94.9", "Russia"),
        entry("95.0", "Russia"),
        entry("95.1", "Russia"),
        entry("95.2", "Russia"),
        entry("95.3", "Russia"),
        entry("95.4", "Russia"),
        entry("95.5", "Poland"),
        entry("95.6", "Poland"),
        entry("95.7", "Slovakia"),
        entry("95.8", "Poland"),
        entry("95.9", "Poland"),
        entry("96.0", "Poland"),
        entry("96.1", "Ukraine"),
        entry("96.2", "Ukraine"),
        entry("96.3", "Ukraine"),
        entry("96.4", "Poland"),
        entry("96.5", "Ukraine"),
        entry("96.6", "Ukraine"),
        entry("96.7", "Poland"),
        entry("96.8", "Ukraine"),
        entry("96.9", "Hungary"),
        entry("97.0", "Serbia"),
        entry("97.1", "Montenegro"),
        entry("97.2", "Montenegro"),
        entry("97.3", "Montenegro"),
        entry("97.4", "Montenegro"),
        entry("97.5", "the Atlantic Ocean"),
        entry("97.6", "the Atlantic Ocean"),
        entry("97.7", "the Atlantic Ocean"),
        entry("97.8", "the Atlantic Ocean"),
        entry("97.9", "the Atlantic Ocean"),
        entry("98.0", "the Atlantic Ocean"),
        entry("98.1", "the Atlantic Ocean"),
        entry("98.2", "the Atlantic Ocean"),
        entry("98.3", "the Atlantic Ocean"),
        entry("98.4", "the Atlantic Ocean"),
        entry("98.5", "the Atlantic Ocean"),
        entry("98.6", "the Atlantic Ocean"),
        entry("98.7", "the Atlantic Ocean"),
        entry("98.8", "the Atlantic Ocean"),
        entry("98.9", "Libya"),
        entry("99.0", "Libya"),
        entry("99.1", "Libya"),
        entry("99.2", "Libya"),
        entry("99.3", "Libya"),
        entry("99.4", "Libya"),
        entry("99.5", "Denmark"),
        entry("99.6", "Egypt"),
        entry("99.7", "the Indian Ocean"),
        entry("99.8", "Saudi Arabia"),
        entry("99.9", "Saudi Arabia"),
        entry("100.0", "Canada")
    );

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
        double lat = -1;
        double lon = -1;
        int radius = -1;
        boolean populationAndRadiusInitialized = false;
        final String foundPercentCirclesFilename = 
            "/home/alexmijo/PopulationCircles/foundPercentageCircles2020.txt";
        // TODO: Modularize into a function
        try(BufferedReader br = new BufferedReader(new FileReader(foundPercentCirclesFilename))) {
            String line = br.readLine();
            while (line != null) {
                String[] foundPercentCircle = line.split(" ");
                if (Math.abs(Double.parseDouble(foundPercentCircle[0]) 
                    - Double.parseDouble(percent)) < 0.0001) {
                    population = (long)Double.parseDouble(foundPercentCircle[4]);
                    lat = Double.parseDouble(foundPercentCircle[3]);
                    lon = Double.parseDouble(foundPercentCircle[2]);
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
        int verySmallFontHeight = 59;
        int veryVerySmallFontHeight = 35;
        double lineSpacing = 1.3;
        int width = old.getWidth();
        int height = old.getHeight() + (int)(bigFontHeight * lineSpacing) 
                     + (int)(veryVerySmallFontHeight * lineSpacing);
        BufferedImage img = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g2d = img.createGraphics();
        g2d.setColor(Color.BLACK);
        g2d.fillRect(0, 0, width, height);
        g2d.setPaint(Color.WHITE);
        g2d.setFont(new Font("Serif", Font.BOLD, bigFontHeight));
        NumberFormat commasFormat = NumberFormat.getInstance();
        commasFormat.setGroupingUsed(true); // this will also round numbers, 3 decimal places
        String title = "The smallest possible circle containing " + percent
                       + "% of the world's population";
        String populationString = commasFormat.format(population) + " people";
        String radiusString = "Radius: " + commasFormat.format(radius) +" km";
        // ImageObserver not needed
        g2d.drawImage(old, 0, height - old.getHeight() 
                      - (int)(veryVerySmallFontHeight * lineSpacing), null);
        g2d.drawString(title, 5, bigFontHeight);
        g2d.setFont(new Font("Serif", Font.BOLD, smallFontHeight));
        g2d.drawString(populationString, 5, smallFontHeight + (int)(bigFontHeight * lineSpacing));
        g2d.drawString(radiusString, 5, smallFontHeight + (int)(smallFontHeight * lineSpacing)
                       + (int)(bigFontHeight * lineSpacing));
        
        g2d.setFont(new Font("Serif", Font.BOLD, verySmallFontHeight));
        int countryY = height - (int)(verySmallFontHeight * (lineSpacing - 1)) - 3;
        String centerString = "Center:";
        DecimalFormat latLonFormatter = new DecimalFormat("#0.000");
        String centerCoords = "(" + latLonFormatter.format(lat) + "°N, " 
                              + latLonFormatter.format(lon) + "°E)";
        if (lon < 0) {
            centerCoords = "(" + latLonFormatter.format(lat) + "°N, " 
                           + latLonFormatter.format(-lon) + "°W)";
        }
        String country = "in " + countries.get(percent);
        g2d.drawString(country, 5, countryY);
        g2d.drawString(centerCoords, 5, countryY - (int)(verySmallFontHeight * lineSpacing));
        g2d.drawString(centerString, 5, countryY - (int)(verySmallFontHeight * lineSpacing) 
                       - (int)(verySmallFontHeight * lineSpacing));
        
        g2d.setFont(new Font("Serif", Font.BOLD, veryVerySmallFontHeight));
        g2d.setPaint(Color.WHITE);
        FontMetrics fm = g2d.getFontMetrics();
        String dataSourceString = 
            "https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-rev11";
        int dataSourceX = width - fm.stringWidth(dataSourceString) - 5;
        int dataSourceY = height - (int)(veryVerySmallFontHeight * (lineSpacing - 1)) - 3;
        g2d.drawString(dataSourceString, dataSourceX, dataSourceY);
        String dataResolutionString = "30 arcsecond resolution";
        int dataResolutionX = width - fm.stringWidth(dataResolutionString) - 5;
        int dataResolutionY = dataSourceY - (int)(veryVerySmallFontHeight * lineSpacing);
        g2d.drawString(dataResolutionString, dataResolutionX, dataResolutionY);
        // g2d.drawString(dataResolutionString, dataResolutionX, dataSourceY);
        String dataYearString = "2020 population data";
        int dataYearX = width - fm.stringWidth(dataYearString) - 5;
        int dataYearY = dataResolutionY - (int)(veryVerySmallFontHeight * lineSpacing);
        g2d.drawString(dataYearString, dataYearX, dataYearY);
        // g2d.drawString(dataYearString, dataYearX, dataResolutionY);
        g2d.dispose();
        return img;
    }

    public static void main(String[] args) {
        for (circleNum = 899; circleNum <= 916; circleNum++) {
            percent = (circleNum / 10) + "." + (circleNum % 10);
            String inputImageFileName =
                "/mnt/c/Users/Administrator/Desktop/linuxPythonMadePercentMaps2020/" + percent 
                + "PercentCircle.png";
            String outputImageFileName =
                "/mnt/c/Users/Administrator/Desktop/linuxPythonMadePercentMaps2020/WithText/0" + circleNum 
                + ".png";
            manipulateImage(inputImageFileName, outputImageFileName);
        }
    }
}