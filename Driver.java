/*
    Author: Rogelio Schevenin Jr.
    Course: CS-310 Data Structures
    Program 4: Genetic Algorithms
    Date: Dec 9, 2020
 */

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class Driver {

    public static Map<String, int[]> matrix = new LinkedHashMap<>();
    public static boolean SHOW_DETAILS;
    public static boolean PAUSE;
    public static final int MUTATION_ODDS = 5;
    public static final int BREED_BY_FITNESS_ODDS = 97;

    public static void main(String[] args) {

        SHOW_DETAILS = args[1].equals("true");
        PAUSE = args[2].equals("true");

        // exit if unable to read file
        if (!readFile(args[0])) {
            System.exit(0);
        }

        // until the algorithm completes, rerun algorithm
        while (!runAlgorithm()) {
            System.out.println("Running algorithm again.");
        }
    }

    // readFile: take file from input and create a matrix/graph out of it
    public static boolean readFile(String file) {
        // Find file from args
        File text = new File(file);

        // Check file
        Scanner scanner;
        try {
            scanner = new Scanner(text);

            String line, location, numbers;
            String[] parts;
            int[] costs;
            int count = 0;

            // Extract from input file
            while (scanner.hasNext()) {
                line = scanner.nextLine();

                // format line
                line = line.replaceAll(", ", " ");
                line = line.replaceAll(",", " ");

                // store location and numbers
                location = line.substring(0, 1);
                numbers = line.substring(2);

                // create array out of numbers
                parts = numbers.split(" ");
                costs = new int[parts.length];

                // put parts into costs as integers
                for (int n = 0; n < parts.length; n++) {
                    costs[n] = Integer.parseInt(parts[n]);
                }

                // add to matrix, create graph
                matrix.put(location, costs);
                System.out.print(location + " -> ");
                System.out.println(Arrays.toString(costs));

                // count amount of locations
                count++;
            }
            scanner.close();
            System.out.println("Successfully loaded adjacency matrix!\nFound " + count + " locations in file.");

            // Displaying information about file to indicate successful graph construction
            System.out.println("Algorithm Test: the cost of going from a to c is " + calculateFitness(matrix, "a", "c") + ".");
            System.out.println("Displaying more information: " + SHOW_DETAILS);
            System.out.println("Pause after each epoch: " + PAUSE + "\n");

            return true;
        } catch (FileNotFoundException e) {
            System.out.println("File not found! Please specify the file's directory (path).");
            return false;
        }
    }

    // runAlgorithm: run the genetic algorithm
    public static boolean runAlgorithm() {
        Scanner scanner = new Scanner(System.in);
        int chromosomes, epochs;

        try {
            System.out.println("Enter a number of initial chromosomes to generate: ");
            // number generated by initial population
            chromosomes = scanner.nextInt();
            System.out.println("Enter a number of epochs (generations) to test: ");
            // number of breeding cycles (crossover) operations
            epochs = scanner.nextInt();

            GeneticAlgorithm algorithm = new GeneticAlgorithm(matrix, chromosomes, epochs, MUTATION_ODDS, BREED_BY_FITNESS_ODDS, SHOW_DETAILS, PAUSE);
            return true;
        } catch (InputMismatchException exception) {
            System.out.println("Could not understand your input. Verify you input is an integer value.\n");
            return false;
        }
    }

    // HELPER METHOD FOR TEST FOLLOWING SUCCESSFUL MATRIX MAPPING
    public static int calculateFitness(Map<String, int[]> matrix, String a, String b) {
        String[] locations = matrix.keySet().toArray(new String[0]);

        // get costs of location a
        int[] costs = matrix.get(a);
        int i = 0, posB = 0;

        // determine the cost of going from a to b
        for (String location : locations) {
            if (location.equals(b)) {
                posB = i;
            }
            i++;
        }

        return costs[posB];
    }
}
