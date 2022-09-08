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


/**
 * Based off of https://stackoverflow.com/questions/2658663
 */
class ImageManipulation {

    private static int percent;

    private static Map<Integer, String> countries = Map.ofEntries(
        entry(1, "Bangladesh"),
        entry(2, "Bangladesh"),
        entry(3, "Bangladesh"),
        entry(4, "India"),
        entry(5, "India"),
        entry(6, "India"),
        entry(7, "India"),
        entry(8, "India"),
        entry(9, "India"),
        entry(10, "India"),
        entry(11, "India"),
        entry(12, "India"),
        entry(13, "India"),
        entry(14, "India"),
        entry(15, "India"),
        entry(16, "India"),
        entry(17, "India"),
        entry(18, "India"),
        entry(19, "India"),
        entry(20, "India"),
        entry(21, "India"),
        entry(22, "India"),
        entry(23, "India"),
        entry(24, "India"),
        entry(25, "China"),
        entry(26, "China"),
        entry(27, "China"),
        entry(28, "China"),
        entry(29, "China"),
        entry(30, "China"),
        entry(31, "China"),
        entry(32, "China"),
        entry(33, "China"),
        entry(34, "China"),
        entry(35, "China"),
        entry(36, "China"),
        entry(37, "India"),
        entry(38, "India"),
        entry(39, "India"),
        entry(40, "India"),
        entry(41, "Myanmar"),
        entry(42, "Myanmar"),
        entry(43, "Myanmar"),
        entry(44, "India"),
        entry(45, "India"),
        entry(46, "Myanmar"),
        entry(47, "Myanmar"),
        entry(48, "India"),
        entry(49, "China"),
        entry(50, "China"),
        entry(51, "Vietnam"),
        entry(52, "China"),
        entry(53, "China"),
        entry(54, "China"),
        entry(55, "China"),
        entry(56, "China"),
        entry(57, "China"),
        entry(58, "Kyrgyzstan"),
        entry(59, "China"),
        entry(60, "Tajikistan"),
        entry(61, "Kyrgyzstan"),
        entry(62, "Tajikistan"),
        entry(63, "Tajikistan"),
        entry(64, "Kazakhstan"),
        entry(65, "Kazakhstan"),
        entry(66, "Uzbekistan"),
        entry(67, "Turkmenistan"),
        entry(68, "Afghanistan"),
        entry(69, "Afghanistan"),
        entry(70, "Turkmenistan"),
        entry(71, "Iran"),
        entry(72, "Iran"),
        entry(73, "Pakistan"),
        entry(74, "Pakistan"),
        entry(75, "Pakistan"),
        entry(76, "Pakistan"),
        entry(77, "Iran"),
        entry(78, "the Indian Ocean"),
        entry(79, "the Indian Ocean"),
        entry(80, "the Indian Ocean"),
        entry(81, "the Indian Ocean"),
        entry(82, "Afghanistan"),
        entry(83, "Pakistan"),
        entry(84, "Pakistan"),
        entry(85, "Iran"),
        entry(86, "Russia"),
        entry(87, "Russia"),
        entry(88, "Russia"),
        entry(89, "Russia"),
        entry(90, "Sweden"),
        entry(91, "Sweden"),
        entry(92, "Germany"),
        entry(93, "Czechia"),
        entry(94, "the Arctic Ocean"),
        entry(95, "Russia"),
        entry(96, "Poland"),
        entry(97, "Serbia"),
        entry(98, "the Atlantic Ocean"),
        entry(99, "Libya"),
        entry(100, "")
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
            "/home/alexmijo/PopulationCircles/foundPercentageCircles.txt";
        // TODO: Modularize into a function
        try(BufferedReader br = new BufferedReader(new FileReader(foundPercentCirclesFilename))) {
            String line = br.readLine();
            while (line != null) {
                String[] foundPercentCircle = line.split(" ");
                if (Integer.parseInt(foundPercentCircle[0]) == percent) {
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
        int verySmallFontHeight = 53;
        int veryVerySmallFontHeight = 35;
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
        int countryY = height - (int)(verySmallFontHeight * (lineSpacing - 1)) - 5;
        String centerString = "Center:";
        DecimalFormat latLonFormatter = new DecimalFormat("#0.000");
        String centerCoords = "(" + latLonFormatter.format(lat) + "°N, " 
                              + latLonFormatter.format(lon) + "°E)";
        String country = "in " + countries.get(percent);
        g2d.drawString(country, 5, countryY);
        g2d.drawString(centerCoords, 5, countryY - (int)(verySmallFontHeight * lineSpacing));
        g2d.drawString(centerString, 5, countryY - (int)(verySmallFontHeight * lineSpacing) 
                       - (int)(verySmallFontHeight * lineSpacing));
        
        g2d.setFont(new Font("Serif", Font.BOLD, veryVerySmallFontHeight));
        g2d.setPaint(Color.WHITE);
        FontMetrics fm = g2d.getFontMetrics();
        String dataSourceString = "https://ghsl.jrc.ec.europa.eu/ghs_pop2019.php";
        int dataSourceX = width - fm.stringWidth(dataSourceString) - 5;
        int dataSourceY = height - (int)(veryVerySmallFontHeight * (lineSpacing - 1)) - 5;
        g2d.drawString(dataSourceString, dataSourceX, dataSourceY);
        String dataResolutionString = "30 arcsecond resolution";
        int dataResolutionX = width - fm.stringWidth(dataResolutionString) - 5;
        int dataResolutionY = dataSourceY - (int)(veryVerySmallFontHeight * lineSpacing);
        g2d.drawString(dataResolutionString, dataResolutionX, dataResolutionY);
        String dataYearString = "2015 population data";
        int dataYearX = width - fm.stringWidth(dataYearString) - 5;
        int dataYearY = dataResolutionY - (int)(veryVerySmallFontHeight * lineSpacing);
        g2d.drawString(dataYearString, dataYearX, dataYearY);
        g2d.dispose();
        return img;
    }

    public static void main(String[] args) {
        for (percent = 1; percent <= 99; percent++) {
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