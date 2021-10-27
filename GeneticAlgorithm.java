/*
    Author: Rogelio Schevenin Jr.
    Course: CS-310 Data Structures
    Program 4: Genetic Algorithms
    Date: Dec 9, 2020
 */

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;

public class GeneticAlgorithm {
    private static boolean SHOW_DETAILS;
    public static int currentEpoch;
    public static String epochSpecifier = "[info-epoch-" + currentEpoch + "]";

    // Constructor
    public GeneticAlgorithm(Map<String, int[]> matrix, int chromosomes, int epochs, boolean output, int mutationOdds, int breedByFitnessOdds) {
        Map<String[], Integer> epochGenome = null;
        String[] bestChromosome = new String[0];
        int bestFitness = 0, totalCost = 0;
        SHOW_DETAILS = output;

        /*
        Done + As the program progresses, it shall write all output to the terminal window (System.out or std::cout <<).
        Done + At the start of each epoch, the program shall display information about the current chromosome population.
        Done + These data shall include the epoch number, the length of the 'most fit' path, and the average fitness of all the chromosomes in the epoch.
        Done + The program may provide a prompt informing the user of the crossover step, but this behavior is optional and may not add any value.
        Done + The program should report when it performs a mutation, or the number of mutations performed, on each epoch.
        Done + This will let you know the operation triggered and verify its frequency.
         */

        // for each generation/epoch
        for (int epoch = 1; epoch <= epochs; epoch++) {
            System.out.println();
            System.out.println("=== Epoch " + epoch + " ===");
            currentEpoch = epoch;

            // if initial epoch, generate initial population
            if (epoch == 1) {
                // all locations in matrix
                LinkedList<String> locations = new LinkedList<>(matrix.keySet());
                // generating random genome ranked by fitness (chromosomes, fitness number)
                Map<String[], Integer> genome = rankFitness(matrix, generateGenome(chromosomes, locations));
                // epoch genome after mutation of crossovers
                epochGenome = mutate(matrix, crossover(matrix, genome, breedByFitnessOdds), mutationOdds);
            } else {
                // epoch genome after mutation of crossovers
                epochGenome = mutate(matrix, crossover(matrix, epochGenome, breedByFitnessOdds), mutationOdds);
            }

            // chromosomes and their fitness from epoch
            LinkedList<String[]> chromosomesThisEpoch = new LinkedList<>(epochGenome.keySet());
            LinkedList<Integer> fitnessesThisEpoch = new LinkedList<>(epochGenome.values());

            // statistics variables
            int chromosomeCount = 0, bestCostThisEpoch = 0, totalCostThisEpoch = 0;
            String[] bestChromosomeThisEpoch = new String[0];


            // for chromosomes in genome
            for (String[] chromosome : chromosomesThisEpoch) {

                // add up all chromosome costs
                totalCostThisEpoch += fitnessesThisEpoch.get(chromosomeCount);
                totalCost += fitnessesThisEpoch.get(chromosomeCount);

                // if the first iteration
                if (chromosomeCount == 0) {
                    bestChromosomeThisEpoch = chromosome;
                    bestCostThisEpoch = fitnessesThisEpoch.get(chromosomeCount);
                }

                // increase count of chromosomes
                chromosomeCount++;
            }

            // if first epoch, initialize the best chromosome and fitness
            // else if current epoch best is better than overall best, replace the overall best
            if (epoch == 1) {
                bestChromosome = bestChromosomeThisEpoch;
                bestFitness = bestCostThisEpoch;
            } else if (bestFitness > bestCostThisEpoch) {
                bestChromosome = bestChromosomeThisEpoch;
                bestFitness = bestCostThisEpoch;
            }

            System.out.println("[info-epoch-" + currentEpoch + "] Most fit path: " + Arrays.toString(bestChromosomeThisEpoch) + " (Cost: " + bestCostThisEpoch + ", Average: " + (totalCostThisEpoch / chromosomeCount) + ")");

            /* pause after each epoch until user input
            if (SHOW_DETAILS) {
                System.out.println("PRESS ENTER TO CONTINUE");
                Scanner s = new Scanner(System.in);
                String input = s.nextLine();
            } else {
                System.out.println();
            }
            */
        }

        if (SHOW_DETAILS) displayChromosomeCost(matrix, bestChromosome);
        System.out.println("[all-epochs] Most fit path: " + Arrays.toString(bestChromosome) + " (Cost: " + bestFitness + ")");

        /*
        Done + The program shall terminate after it completes the number of steps established by the user when the program launched.
        Done + Upon exit, it shall display information about the winning chromosome (path).
        Done + It shall include the sequence of nodes to visit as well as the path's length.
         */
    }

    /*
    Done + Generate an Initial Population: Randomly generate several complete sequences of nodes to visit.
    Done + If there are twenty-five nodes, the sequence shall include each one exactly once.
    Done + Keeping with the 'genetic' analogy, these sequences are the chromosomes for the algorithm.
     */
    public static LinkedList<String[]> generateGenome(int n, LinkedList<String> locations) {
        System.out.println(epochSpecifier + " Generating initial genome.");
        LinkedList<String[]> chromosomes = new LinkedList<>();

        // make n chromosomes
        for (int i = 0; i < n; i++) {
            String[] chromosome;

            // add first random chromosome to linked list (nothing to compare to)
            if (i == 0) {
                chromosome = randomChromosome(locations);
                chromosomes.add(chromosome);

                // output
                if (SHOW_DETAILS) System.out.println(epochSpecifier + " Generated random chromosome -> " + Arrays.toString(chromosome) + " (" + i + "/" + n + ")");
            } else {
                // create new random chromosome
                chromosome = randomChromosome(locations);

                // ensure it's not already in list
                while (containsChromosome(chromosome, chromosomes)) {
                    chromosome = randomChromosome(locations);
                }

                // add to list
                chromosomes.add(chromosome);

                //output
                if (SHOW_DETAILS) System.out.println(epochSpecifier + " Generated random chromosome -> " + Arrays.toString(chromosome) + " (" + i + "/" + n + ")");
            }
        }

        // return chromosomes (genome)
        return chromosomes;
    }

    /*
    Done + Fitness Function: Survival of the fittest requires us to identify a metric with which we can rank each chromosome.
    Done + For this problem, the fitness of a solution should reflect the length of the chromosome's path.
    Done + For us, smaller numbers are more fit.
     */
    public static Map<String[], Integer> rankFitness(Map<String, int[]> matrix, LinkedList<String[]> chromosomes) {
        Map<String[], Integer> rankedFitness = new LinkedHashMap<>();

        System.out.println(epochSpecifier + " Analyzing fitness of chromosomes.");

        // traverse all chromosomes
        for (String[] chromosome : chromosomes) {
            int totalFitness = 0;

            // traverse all locations in chromosome
            for (int i = 0; i < chromosome.length; i++) {
                // if not last in chromosome
                if (i != chromosome.length - 1) {
                    totalFitness += calculateFitness(matrix, chromosome[i], chromosome[i + 1]);
                } else {
                    totalFitness += calculateFitness(matrix, chromosome[i], chromosome[i]);
                }
            }

            // add last route (from last location to first location) to total fitness of chromosome
            totalFitness += calculateFitness(matrix, chromosome[chromosome.length - 1], chromosome[0]);

            if (SHOW_DETAILS) System.out.println(epochSpecifier + " Analyzed chromosome fitness -> " + Arrays.toString(chromosome) + ": " + totalFitness);

            // add sequence and the cost of the sequence to the map
            rankedFitness.put(chromosome, totalFitness);
        }

        // return sorted by best fitness (lowest number) to worst fitness (highest number)
        return sortByFitness(rankedFitness);
    }

    /*
    Done + Crossover: Based on their fitness, partner chromosomes up for breeding.
    Done + The strategy for partnering up chromosomes varies with each solution, but generally the more fit chromosomes have a greater chance of breeding than the less fit ones, but research suggests including a few of the less-fit chromosomes in the breeding cycle can lead to a better solution.
     */
    public static Map<String[], Integer> crossover(Map<String, int[]> matrix, Map<String[], Integer> genome, int breedByFitnessOdds) {
        LinkedList<String[]> chromosomes = new LinkedList<>(genome.keySet());

        System.out.println(epochSpecifier + " Crossover: Partnering chromosomes for breeding.");

        // divide best half of chromosomes into fit list
        LinkedList<String[]> fitChromosomes = new LinkedList<>();
        for (int i = 0; i < chromosomes.size() / 2; i++) {
            fitChromosomes.add(chromosomes.get(i));
            if (SHOW_DETAILS) System.out.println(epochSpecifier + " Chromosome classified as fit -> " + Arrays.toString(chromosomes.get(i)));
        }

        System.out.println(epochSpecifier + " Crossover: Successfully generated fit chromosomes list.");

        // divide best half of chromosomes into non-fit list
        LinkedList<String[]> nonFitChromosomes = new LinkedList<>();
        for (int i = chromosomes.size() / 2; i < chromosomes.size(); i++) {
            nonFitChromosomes.add(chromosomes.get(i));
            if (SHOW_DETAILS) System.out.println(epochSpecifier + " Chromosome classified as non-fit -> " + Arrays.toString(chromosomes.get(i)));
        }

        System.out.println(epochSpecifier + " Crossover: Successfully generated non-fit chromosomes list.");

        System.out.println(epochSpecifier + " Crossover: Breeding chromosomes.");

        // list for chromosomes after preforming crossovers
        LinkedList<String[]> newChromosomes = new LinkedList<>();

        // until newChromosome list has been filled with chromosomes
        while (newChromosomes.size() != chromosomes.size()) {

            Random r = new Random();

            String[] firstParent;
            String[] secondParent;
            String[] c;

            String type;

            // top 10% of fit chromosomes, unless initial population is too small for top 10%, then random from fit list
            int best = (int) (((int) (fitChromosomes.size() * (0.1)) == 0) ? (fitChromosomes.size()) : ((fitChromosomes.size() * (0.1))));

            // random chance that unique non-fit parents breed
            int chance = r.nextInt(100);

            if (chance < breedByFitnessOdds) {
                do {
                    firstParent = fitChromosomes.get(r.nextInt(best));
                    secondParent = fitChromosomes.get(r.nextInt(best));
                } while (Arrays.equals(firstParent, secondParent));
                type = "fit";

                if (SHOW_DETAILS) System.out.println(epochSpecifier + " Parent 1 (" + type + ") -> " + Arrays.toString(firstParent) + " + Parent 2 (" + type + ") -> " + Arrays.toString(secondParent));
            } else {
                do {
                    firstParent = nonFitChromosomes.get(r.nextInt(nonFitChromosomes.size()));
                    secondParent = nonFitChromosomes.get(r.nextInt(nonFitChromosomes.size()));
                } while (Arrays.equals(firstParent, secondParent));
                type = "non-fit";
            }

            c = crossoverChromosomes(firstParent, secondParent);

            if (SHOW_DETAILS) System.out.println(epochSpecifier + " Bred child chromosome -> " + Arrays.toString(c) + " (" + newChromosomes.size() + "/" + chromosomes.size() + ")");

            if (!containsChromosome(c, newChromosomes)) {
                newChromosomes.add(c);
            }
        }

        System.out.println(epochSpecifier + " Successfully performed crossovers.");

        // rank them by fitness and return the genome
        return rankFitness(matrix, newChromosomes);
    }

    /*
    Done + Mutation: Each resultant chromosome (the output of the crossover process) has a chance to randomly mutate.
    Done + For a mutation, simply flip the positions of two nodes in the sequence.
     */
    public static Map<String[], Integer> mutate(Map<String, int[]> matrix, Map<String[], Integer> genome, int mutationOdds) {
        LinkedList<String[]> chromosomes = new LinkedList<>(genome.keySet());
        LinkedList<String[]> newChromosomes = new LinkedList<>();

        // adjust mutation odds for smaller variance
        if (chromosomes.size() <= 10) {
            mutationOdds *= 15;
        } else if (chromosomes.size() <= 100) {
            mutationOdds *= 5;
        }

        // go through old list of chromosomes and put them through mutation process
        for (String[] chromosome : chromosomes) {
            newChromosomes.add(mutateChromosome(chromosome, mutationOdds));
        }

        // ensure list of chromosomes did not change
        if (chromosomes.size() == newChromosomes.size()) {
            System.out.println(epochSpecifier + " Successfully performed mutations.");
            return rankFitness(matrix, newChromosomes);
        } else {
            System.out.println(epochSpecifier + " Mutation error: size is different for old chromosomes and new chromosomes!");
            return null;
        }
    }

    /*
    HELPER METHODS
     */

    // sortByFitness: sorts genome (ascending order) based on fitness
    private static Map<String[], Integer> sortByFitness(Map<String[], Integer> genome) {
        // collections sort list by value
        List<Map.Entry<String[], Integer>> list = new ArrayList<>(genome.entrySet());
        list.sort(Map.Entry.comparingByValue());

        // put back into map
        Map<String[], Integer> sorted = new LinkedHashMap<>();
        for (Map.Entry<String[], Integer> entry : list) {
            sorted.put(entry.getKey(), entry.getValue());
        }

        return sorted;
    }

    // randomChromosome: generates a random chromosome
    private static String[] randomChromosome(LinkedList<String> locations) {
        Random r = new Random();

        // array for a chromosome with enough space to fit exactly one of every location
        String[] chromosome = new String[locations.size()];

        // based on amount of locations
        for (int x = 0; x < locations.size(); x++) {
            int random = r.nextInt(locations.size());

            // while chromosome already contains random location
            while (Arrays.asList(chromosome).contains(locations.get(random))) {
                // generate new location to add next to chromosome
                random = r.nextInt(locations.size());
            }

            // add random location to sequence
            chromosome[x] = (locations.get(random));
        }

        // return chromosome
        return chromosome;
    }

    // crossoverChromosomes: performs a single crossover on two parent chromosomes
    private static String[] crossoverChromosomes(String[] a, String[] b) {
        Random r = new Random();

        String[] parentA = new String[b.length];
        String[] parentB = new String[b.length];
        String[] child = new String[b.length];
        String[] c = new String[child.length];

        // clone chromosome
        for (int i = 0; i < parentA.length; i++) {
            parentA[i] = a[i];
            parentB[i] = b[i];
            child[i] = parentB[i];
        }

        // do size/2 swaps to create child of parentA and parentB
        for (int x = 0; x < parentA.length / 2; x++) {
            // pick random index to copy
            int swapIndex = r.nextInt(parentA.length);
            int secondSwapIndex = 0;

            // swap values
            String parentValueToSwap = parentA[swapIndex];
            String childValueToSwap = child[swapIndex];

            // find second swap index
            for (int y = 0; y < parentA.length; y++) {
                if (child[y].equals(parentValueToSwap)) {
                    secondSwapIndex = y;
                }
            }

            // perform swap
            child[swapIndex] = parentValueToSwap;
            child[secondSwapIndex] = childValueToSwap;

            // get string value of ints in array
            for (int z = 0; z < child.length; z++) {
                c[z] = String.valueOf(child[z]);
            }
        }

        return c;
    }

    // mutateChromosome: performs a mutation on 2% of chromosomes
    private static String[] mutateChromosome(String[] chromosome, int mutationOdds) {
        Random r = new Random();

        String[] reference = new String[chromosome.length];
        String[] original = new String[chromosome.length];

        // clone chromosome
        for (int i = 0; i < chromosome.length; i++) {
            reference[i] = (chromosome[i]);
            original[i] = reference[i];
        }

        // chance the chromosome mutates
        if (r.nextInt(100) < mutationOdds) {

            // pick random index to copy
            int swapIndex = r.nextInt(chromosome.length);
            int secondSwapIndex = 0;

            // swap values
            String referenceValueToSwap = reference[swapIndex];
            String originalValueToSwap = original[swapIndex];

            // find second swap index
            for (int y = 0; y < chromosome.length; y++) {
                if (original[y].equals(referenceValueToSwap)) {
                    secondSwapIndex = y;
                }
            }

            // perform swap
            original[swapIndex] = referenceValueToSwap;
            original[secondSwapIndex] = originalValueToSwap;

            // get string value of ints in array
            for (int z = 0; z < chromosome.length; z++) {
                chromosome[z] = String.valueOf(original[z]);
            }

            if (SHOW_DETAILS) System.out.println(epochSpecifier + " Mutation: Chromosome mutated -> " + Arrays.toString(chromosome));

        }

        return chromosome;
    }

    // containsChromosome: returns true if a list of chromosomes contains a chromosome
    private static boolean containsChromosome(String[] target, LinkedList<String[]> chromosomes) {

        // sort
        chromosomes.sort(Comparator.comparing(o -> o[1]));

        // search
        for (String[] chromosome : chromosomes) {
            if (Arrays.equals(chromosome, target)) {
                return true;
            }
        }

        return false;
    }

    // calculateFitness: determines fitness of from one location to another
    private static int calculateFitness(Map<String, int[]> matrix, String a, String b) {
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

    // displayCost: prints a breakdown of the cost of the path a chromosome takes
    public static void displayChromosomeCost(Map<String, int[]> matrix, String[] chromosome) {
        int fitness = 0;
        int totalFitness = 0;
        String a = "";
        String b = "";
        int cost = 0;

        // traverse all locations in chromosome
        for (int i = 0; i < chromosome.length; i++) {
            // if not last in chromosome
            if (i != chromosome.length - 1) {
                a = chromosome[i];
                b = chromosome[i+1];
            } else {
                a = chromosome[i];
                b = chromosome[i];
            }

            fitness = calculateFitness(matrix, a, b);
            totalFitness += fitness;
            System.out.println("[" + a + "] to [" + b + "]: $" + calculateFitness(matrix,a,b));
        }

        // add last route (from last location to first location) to total fitness of chromosome
        totalFitness += calculateFitness(matrix, chromosome[chromosome.length - 1], chromosome[0]);
        System.out.println("[" + chromosome[chromosome.length - 1] + "] to [" + chromosome[0] + "]: $" + calculateFitness(matrix, chromosome[chromosome.length - 1], chromosome[0]));
    }
}
