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

import java.util.List;
import java.util.Map;
import java.lang.Math;


/**
 * Based off of https://stackoverflow.com/questions/2658663
 */
class ImageManipulation {

    private static String percent;
    private static int circleNum;

    // May be too big to compile. I split it into thirds and did each third separately. Should put
    //  this in a text file instead in the future (TODO).
    private static List<String> countries = List.of(
    "Palikir, Federated States of Micronesia", 
    "Victoria, Seychelles",
    "Yaren District, Nauru",
    "South Tarawa, Kiribati",
    "Delap-Uliga-Djarrit, Marshall Islands",
    "Funafuti, Tuvalu",
    "Reykjavik, Iceland",
    "Apia, Samoa",
    "Port Vila, Vanuatu",
    "Suva, Fiji",
    "Honiara, Solomon Islands",
    "Nuku'alofa, Tonga",
    "Wellington, New Zealand",
    "Port Louis, Mauritius",
    "Ngerulmud, Palau",
    "Paramaribo, Suriname",
    "Port Moresby, Papua New Guinea",
    "Windhoek, Namibia",
    "Georgetown, Guyana",
    "Canberra, Australia",
    "Praia, Cape Verde",
    "Ulaanbaatar, Mongolia",
    "Sucre, Bolivia",
    "Bridgetown, Barbados",
    "Dili, East Timor",
    "Castries, Saint Lucia",
    "St. John's, Antigua and Barbuda",
    "Kingstown, Saint Vincent and the Grenadines",
    "Roseau, Dominica",
    "Astana, Kazakhstan",
    "Port of Spain, Trinidad and Tobago",
    "St. George's, Grenada",
    "Santiago, Chile",
    "Lima, Peru",
    "Bandar Seri Begawan, Brunei",
    "Antananarivo, Madagascar",
    "Basseterre, Saint Kitts and Nevis",
    "Nouakchott, Mauritania",
    "Havana, Cuba",
    "Dakar, Senegal",
    "Port-au-Prince, Haiti",
    "San Jose, Costa Rica",
    "Montevideo, Uruguay",
    "Luanda, Angola",
    "Buenos Aires, Argentina",
    "Nassau, Bahamas",
    "Kingston, Jamaica",
    "Banjul, Gambia",
    "Muscat, Oman",
    "Managua, Nicaragua",
    "Santo Domingo, Dominican Republic",
    "Tegucigalpa, Honduras",
    "Bissau, Guinea-Bissau",
    "Tripoli, Libya",
    "Quito, Ecuador",
    "Caracas, Venezuela",
    "Abu Dhabi, United Arab Emirates",
    "San Salvador, El Salvador",
    "Asuncion, Paraguay",
    "Maseru, Lesotho",
    "Moroni, Comoros",
    "Mogadishu, Somalia",
    "Bangui, Central African Republic",
    "Kinshasa, Democratic Republic of the Congo",
    "Panama City, Panama",
    "Belmopan, Belize",
    "Brazzaville, Republic of the Congo",
    "Guatemala City, Guatemala",
    "Doha, Qatar",
    "Lobamba, Eswatini",
    "Maputo, Mozambique",
    "Gaborone, Botswana",
    "Conakry, Guinea",
    "Freetown, Sierra Leone",
    "Pretoria, South Africa",
    "Monrovia, Liberia",
    "Bogota, Colombia",
    "Valletta, Malta",
    "Riyadh, Saudi Arabia",
    "Lisbon, Portugal",
    "Manama, Bahrain",
    "Helsinki, Finland",
    "Rabat, Morocco",
    "Lusaka, Zambia",
    "Oslo, Norway",
    "Stockholm, Sweden",
    "Tallinn, Estonia",
    "Male, Maldives",
    "Bishkek, Kyrgyzstan",
    "Manila, Philippines",
    "Ashgabat, Turkmenistan",
    "Khartoum, Sudan",
    "Brasilia, Brazil",
    "Ottawa, Canada",
    "Lilongwe, Malawi",
    "Sanaa, Yemen",
    "Harare, Zimbabwe",
    "Moscow, Russia",
    "Kuala Lumpur, Malaysia",
    "Mexico City, Mexico",
    "Bamako, Mali",
    "Tokyo, Japan",
    "Tunis, Tunisia",
    "Baku, Azerbaijan",
    "Athens, Greece",
    "Algiers, Algeria",
    "Kuwait City, Kuwait",
    "Madrid, Spain",
    "N’Djamena, Chad".
    "Libreville, Gabon",
    "Dublin, Ireland",
    "Yamoussoukro, Ivory Coast",
    "Tashkent, Uzbekistan",
    "Tehran, Iran",
    "Tbilisi, Georgia",
    "Sao Tome, Sao Tome and Principe",
    "Riga, Latvia",
    "Singapore",
    "Washington, D.C., USA",
    "Djibouti City, Djibouti",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
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
        String title = percent
                       + "th: " + ;
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
        for (circleNum = 1; circleNum <= 1000; circleNum++) {
            percent = (circleNum / 10) + "." + (circleNum % 10);
            String inputImageFileName =
                "/mnt/d/PopCirclesStuff/reprojected/pngs/" + circleNum 
                + ".png";
            String outputImageFileName;
            if (circleNum < 10) {
                outputImageFileName =
                    "/mnt/d/PopCirclesStuff/WithText/00" + circleNum 
                    + ".png";
            } else if (circleNum < 100) {
                outputImageFileName =
                    "/mnt/d/PopCirclesStuff/WithText/0" + circleNum 
                    + ".png";
            } else {
                outputImageFileName =
                    "/mnt/d/PopCirclesStuff/WithText/" + circleNum 
                    + ".png";
            }
            manipulateImage(inputImageFileName, outputImageFileName);
        }
    }
}