import cern.jet.random.Binomial;
import cern.jet.random.Exponential;
import cern.jet.random.Uniform;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;


public class stochasticSIR {

    infections infectionHistory = new infections(20000);
    List<Integer> currentInfected = new ArrayList();
    List<Integer> recovered = new ArrayList();
    List<Integer> exposed = new ArrayList();

    private BufferedWriter writer1 = null;
    private BufferedWriter writer2 = null;

    // simulation parameters

    double maxTime = 360;
    double tau = 1.0;
    int N = 100000;
    int Y_init = 5;
    double incubationPeriod = 4.0;
    double infectiousPeriod = 14; // days
    double R0 = 5;
    double gamma = 1.0 / ((double) infectiousPeriod);
    double delta = 1.0 / ((double) incubationPeriod);

    double beta = R0 / ((double) infectiousPeriod);

    int S_init = (int) Math.floor(N / R0);
    int I_init = (int) Math.floor(N - S_init);
    int R_init = N - S_init - I_init;

    double sample_rate = 5;  // 5 days
    int pool_size = 100; // 50 individuals
    int n_pools = 3;

    DescriptiveStatistics ctStats = new DescriptiveStatistics();
    DescriptiveStatistics readStats = new DescriptiveStatistics();

    public void runSEIR() {

        File outputfile = new File("infections_through_time.csv");

        try {
            writer1 = new BufferedWriter(new FileWriter(outputfile));
        } catch (IOException e) {
            e.printStackTrace();
        }


        int Y_curr = 0;
        int E_curr = 5;
        int X_curr = N - E_curr; // initially assuming no recovered individuals
        int Z_curr = N - X_curr - E_curr - Y_curr; //should be zero at the start of epidemic

        System.out.println("*** Initial States (X, Y, Z) ***");
        System.out.println(X_curr + "\t" + E_curr + "\t" + Y_curr + "\t" + Z_curr);
        System.out.println();

        double t_curr = 0;
        double t_next;


        double BSI = beta * (X_curr / N) * Y_curr;
        double EI = delta * E_curr;
        double gammaI = gamma * Y_curr;


        int seed = (int) System.currentTimeMillis();

        double p_SI = 1 - Math.exp(-beta * (Y_curr / (double) N));
        double p_EI = 1 - Math.exp(-delta);
        double p_IR = 1 - Math.exp(-gamma);

        RandomEngine engine = new MersenneTwister(seed);
        Binomial binomial = new Binomial(10, 0.2, engine);

        initialiseHistory(E_curr, infectionHistory, exposed);

        int current_id = Y_curr + 1;

        t_curr = t_curr + tau;


        while (t_curr < maxTime) {

//            System.out.println(t_curr + "\t" +  X_curr + "\t" + Y_curr + "\t" + E_curr + "\t" + Z_curr);


            int new_E = 0;
            if (p_SI > 0 && X_curr > 0) {
                new_E = binomial.nextInt(X_curr, p_SI);
            }

            int new_I = 0;
            if (p_EI > 0 && E_curr > 0) {
                new_I = binomial.nextInt(E_curr, p_EI);
            }

            int new_Z = 0;
            if (p_IR > 0 && Y_curr > 0) {
                new_Z = binomial.nextInt(Y_curr, p_IR);
            }

            X_curr = X_curr - new_E;
            E_curr = E_curr + new_E - new_I;
            Y_curr = Y_curr + new_I - new_Z;
            Z_curr = Z_curr + new_Z;


            // new exposeds
            for (int i = 0; i < new_E; i++) {

                int new_id = current_id++;
                infectionHistory.logInfection(new_id, t_curr, -100);
                exposed.add(new_id);
            }

            // new infected
            for (int i = 0; i < new_I; i++) {

                int index = Uniform.staticNextIntFromTo(0, exposed.size() - 1);
                Integer infectious = exposed.get(index);
                currentInfected.add(infectious);
                exposed.remove(infectious);

//                infectionHistory.setBirth(index, t_curr);

            }

            // new recovereds
            List<Double> weights = new ArrayList();

//            for(int i = 0; i < new_Z; i++) {
//
//                int index = Uniform.staticNextIntFromTo(0, currentInfected.size()-1);
//                Integer recoveredInfection = currentInfected.get(index);
//                recovered.add(recoveredInfection);
//                currentInfected.remove(recoveredInfection);
//                infectionHistory.setDeath(recoveredInfection, t_curr);
//            }


            for (Integer s : currentInfected) {

                int index = infectionHistory.id.indexOf(s);
                double time_since_infection = t_curr - infectionHistory.getBirth(index);
                weights.add(time_since_infection);
            }


            List<Integer> recoveredInfections = utils.sampling.multinomSampDouble(weights, currentInfected, new_Z);

            for (Integer r : recoveredInfections) {

                infectionHistory.setDeath(r, t_curr);
                currentInfected.remove(r);
                recovered.add(r);

            }

            if (exposed.size() != E_curr) {
                System.out.println("not right " + E_curr + "\t" + exposed.size());
            }


            // update probabilities

            p_SI = 1 - Math.exp(-beta * (Y_curr / (double) N));
            p_IR = 1 - Math.exp(-gamma);


            t_curr += tau;

            if (E_curr == 0) {
                break;
            }

            if (t_curr % 5 == 0) {

                int[] pool = sampleIndividuals(X_curr, E_curr, Y_curr, Z_curr, 8);

                List<Double> CtDist = getCtDistribution(pool, t_curr);

                ctStats.clear();
                readStats.clear();
                CtDist.forEach(Ct -> {
//                    if(Ct < 37.0) {
//                        Ct = Ct + poolDilutionFactor;
//                    }
//                    else{
//                        Ct = 40.0;
//                    }

                    double read_count = (Math.pow(2, -0.47 * Ct + 14.4));
                    if (read_count < 0) {
                        read_count = 0;
                    }
                    ctStats.addValue(Ct);
                    readStats.addValue(Math.round(read_count));
                });

                double mean = Math.round(ctStats.getMean() * 1e4) / 1e4;
                double var = Math.round(ctStats.getStandardDeviation() * 1e4) / 1e4;
                double skew = Math.round(ctStats.getSkewness() * 1e4) / 1e4;
                double pooled_reads = Math.round(readStats.getSum() * 1e4) / 1e4;

                System.out.println(t_curr + "\t" + (pooled_reads / (double) pool_size) + "\t" + skew + "\t" + new_I / (double) N + "\t" + (Y_curr) / (double) N +
                        "\t" + pool[0] + "\t" + pool[1] + "\t" + pool[2]);


            }


        }

    }

    public void runSIR() {

        File outputfile = new File("infections_through_time.csv");

        try {
            writer1 = new BufferedWriter(new FileWriter(outputfile));
        } catch (IOException e) {
            e.printStackTrace();
        }


        int Y_curr = Y_init;
        int X_curr = N - Y_curr; // initially assuming no recovered individuals
        int Z_curr = N - X_curr - Y_curr; //should be zero at the start of epidemic

        System.out.println("*** Initial States (X, Y, Z) ***");
        System.out.println(X_curr + "\t" + Y_curr + "\t" + Z_curr);
        System.out.println();

        double t_curr = 0;
        double t_next;


        double BSI = beta * (X_curr / N) * Y_curr;
        double gammaI = gamma * Y_curr;


        int seed = (int) System.currentTimeMillis();

        double p_SI = 1 - Math.exp(-beta * (Y_curr / (double) N));
        double p_IR = 1 - Math.exp(-gamma);

        RandomEngine engine = new MersenneTwister(seed);
        Binomial binomial = new Binomial(10, 0.2, engine);

        initialiseHistory(Y_curr, infectionHistory, currentInfected);

        int current_id = Y_curr + 1;

        t_curr = t_curr + tau;


        while (t_curr < maxTime) {

//            System.out.println(t_curr + "\t" +  X_curr + "\t" + Y_curr + "\t" + E_curr + "\t" + Z_curr);

            int new_I = binomial.nextInt(X_curr, p_SI);
            int new_Z = binomial.nextInt(Y_curr, p_IR);


            X_curr = X_curr - new_I;
            Y_curr = Y_curr + new_I - new_Z;
            Z_curr = Z_curr + new_Z;


            // new infected
            for (int i = 0; i < new_I; i++) {

                int new_id = current_id++;
                infectionHistory.logInfection(new_id, t_curr, -100);
                currentInfected.add(new_id);
            }


            // new recovereds
            List<Double> weights = new ArrayList();

            for (Integer s : currentInfected) {

                int index = infectionHistory.id.indexOf(s);
                double time_since_infection = t_curr - infectionHistory.getBirth(index);
                weights.add(time_since_infection);
            }


            List<Integer> recoveredInfections = utils.sampling.multinomSampDouble(weights, currentInfected, new_Z);

            for (Integer r : recoveredInfections) {

                int index = infectionHistory.id.indexOf(r);
                infectionHistory.setDeath(index, t_curr);
                currentInfected.remove(r);
                recovered.add(r);

            }

            if (currentInfected.size() != Y_curr) {
                System.out.println("not right " + Y_curr + "\t" + currentInfected.size());
            }


            // update probabilities

            p_SI = 1 - Math.exp(-beta * (Y_curr / (double) N));
            p_IR = 1 - Math.exp(-gamma);


            t_curr += tau;

            if (Y_curr == 0) {
                break;
            }

            if (t_curr % 5 == 0) {

                int[] pool = sampleIndividuals(X_curr, Y_curr, Z_curr, pool_size);

                List<Double> CtDist = getCtDistribution(pool, t_curr);

                ctStats.clear();
                readStats.clear();
                CtDist.forEach(Ct -> {
//                    if(Ct < 37.0) {
//                        Ct = Ct + poolDilutionFactor;
//                    }
//                    else{
//                        Ct = 40.0;
//                    }
                    ctStats.addValue(Ct);

                    double reads = Math.round(Math.pow(2, -0.47 * Ct + 14.4));
                    if (reads < 0) {
                        reads = 0;
                    }
                    readStats.addValue(reads);

                });

                double mean = Math.round(ctStats.getMean() * 1e4) / 1e4;
                double var = Math.round(ctStats.getStandardDeviation() * 1e4) / 1e4;
                double skew = Math.round(ctStats.getSkewness() * 1e4) / 1e4;
                double pooled_reads = Math.round(readStats.getSum() * 1e4) / 1e4;


                System.out.println(t_curr + "\t" + pooled_reads / (double) pool_size + "\t" + skew + "\t" + new_I / (double) N + "\t" + (Y_curr) / (double) N
                        + "\t" + pool[0] + "\t" + pool[1] + "\t" + pool[2]);


            }


        }

    }

    public void runSIR(int n_sims) {

        DescriptiveStatistics ctStats_25 = new DescriptiveStatistics();
        DescriptiveStatistics readStats_25 = new DescriptiveStatistics();
        DescriptiveStatistics ctStats_50 = new DescriptiveStatistics();
        DescriptiveStatistics readStats_50 = new DescriptiveStatistics();
        DescriptiveStatistics ctStats_100 = new DescriptiveStatistics();
        DescriptiveStatistics readStats_100 = new DescriptiveStatistics();

        File outputfile1 = new File("infections_through_time_N_"+N+".csv");

        try {
            writer1 = new BufferedWriter(new FileWriter(outputfile1));
            writer1.write("Time" + "," +
                    "Pooled reads (25)" + "," + "Skew (25) " + "," +
                    "Pooled reads (50)" + "," + "Skew (50) " + "," +
                    "Pooled reads (100)" + "," + "Skew (100) " + "," +
                    "Incidence (pop)" + "," + "Prevalence (pop)" + "," +
                    "S (pop)" + "," + "I (pop)" + "," + "R (pop)" + "," +
                    "S (25)" + "," + "I (25)" + "," + "R (25)" + "," +
                    "S (50)" + "," + "I (50)" + "," + "R (50)" + "," +
                    "S (100)" + "," + "I (100)" + "," + "R (100)" + "," +
                    "pool Ct (25)" + "," + "pool Ct (50)" + "," + "pool Ct (100)" + "," +

                    "Sim no\n");

        } catch (IOException e) {
            e.printStackTrace();
        }


        int Y_curr = Y_init;
        int X_curr = N - Y_curr; // initially assuming no recovered individuals
        int Z_curr = N - X_curr - Y_curr; //should be zero at the start of epidemic

        System.out.println("*** Initial States (X, Y, Z) ***");
        System.out.println(X_curr + "\t" + Y_curr + "\t" + Z_curr);
        System.out.println();

        double t_curr = 0;


        int seed = (int) System.currentTimeMillis();

        double p_SI = 1 - Math.exp(-beta * (Y_curr / (double) N));
        double p_IR = 1 - Math.exp(-gamma);

        RandomEngine engine = new MersenneTwister(seed);
        Binomial binomial = new Binomial(10, 0.2, engine);

        initialiseHistory(Y_curr, infectionHistory, currentInfected);

        int current_id = Y_curr + 1;

        t_curr = t_curr + tau;


        while (t_curr < maxTime) {

//            System.out.println(t_curr + "\t" +  X_curr + "\t" + Y_curr + "\t" + E_curr + "\t" + Z_curr);


            int new_I = binomial.nextInt(X_curr, p_SI);
            int new_Z = binomial.nextInt(Y_curr, p_IR);


            X_curr = X_curr - new_I;
            Y_curr = Y_curr + new_I - new_Z;
            Z_curr = Z_curr + new_Z;


            // new infected
            for (int i = 0; i < new_I; i++) {

                int new_id = current_id++;
                infectionHistory.logInfection(new_id, t_curr, -100);
                currentInfected.add(new_id);
            }


            // new recovereds
            List<Double> weights = new ArrayList();

            for (Integer s : currentInfected) {

                int index = infectionHistory.id.indexOf(s);
                double time_since_infection = t_curr - infectionHistory.getBirth(index);
                weights.add(time_since_infection);
            }


            List<Integer> recoveredInfections = utils.sampling.multinomSampDouble(weights, currentInfected, new_Z);

            for (Integer r : recoveredInfections) {

                int index = infectionHistory.id.indexOf(r);
                infectionHistory.setDeath(index, t_curr);
                currentInfected.remove(r);
                recovered.add(r);

            }

            if (currentInfected.size() != Y_curr) {
                System.out.println("not right " + Y_curr + "\t" + currentInfected.size());
            }


            // update probabilities

            p_SI = 1 - Math.exp(-beta * (Y_curr / (double) N));
            p_IR = 1 - Math.exp(-gamma);


            t_curr += tau;

            if (Y_curr == 0) {
                break;
            }

            if (t_curr % 20 == 0) {

                for (int ns = 0; ns < n_sims; ns++) {

                    int[] pool_25 = sampleIndividuals(X_curr, Y_curr, Z_curr, 25);
                    int[] pool_50 = sampleIndividuals(X_curr, Y_curr, Z_curr, 50);
                    int[] pool_100 = sampleIndividuals(X_curr, Y_curr, Z_curr, 100);

                    List<Double> CtDist_25 = getCtDistribution(pool_25, t_curr);
                    List<Double> CtDist_50 = getCtDistribution(pool_50, t_curr);
                    List<Double> CtDist_100 = getCtDistribution(pool_100, t_curr);

                    ctStats.clear();
                    readStats.clear();


                    CtDist_25.forEach(Ct -> {

                        ctStats.addValue(Ct);

                        double reads = Math.round(Math.pow(2, -1.0 * Ct + 37.0));
                        if (reads < 0) {
                            reads = 0;
                        }
                        double dil_reads = Math.floor(reads/25.0);
                        readStats.addValue(dil_reads);

                    });

                    double mean_25 = Math.round(readStats.getMean() * 25.0 * 1e4) / 1e4;
                    double var_25 = Math.round(readStats.getStandardDeviation() * 1e4) / 1e4;
                    double skew_25 = Math.round(readStats.getSkewness() * 1e4) / 1e4;
                    double pooled_reads_25 = Math.round(readStats.getSum() * 1e4) / 1e4;
                    double pool_Ct_25 = ((Math.log(readStats.getSum())/Math.log(2))-37)*(-1.0);


                    ctStats.clear();
                    readStats.clear();


                    CtDist_50.forEach(Ct -> {

                        ctStats.addValue(Ct);

                        double reads = Math.round(Math.pow(2, -1.0 * Ct + 37.0));
                        if (reads < 0) {
                            reads = 0;
                        }
                        double dil_reads = Math.floor(reads/50.0);
                        readStats.addValue(dil_reads);

                    });


                    double mean_50 = Math.round(readStats.getMean() * 1e4) / 1e4;
                    double var_50 = Math.round(readStats.getStandardDeviation() * 1e4) / 1e4;
                    double skew_50 = Math.round(readStats.getSkewness() * 1e4) / 1e4;
                    double pooled_reads_50 = Math.round(readStats.getSum() * 1e4) / 1e4;
                    double pool_Ct_50 = ((Math.log(readStats.getSum())/Math.log(2))-37)*(-1.0);

                    ctStats.clear();
                    readStats.clear();


                    CtDist_100.forEach(Ct -> {

                        ctStats.addValue(Ct);

                        double reads = Math.round(Math.pow(2, -1.0 * Ct + 37.0));
                        if (reads < 0) {
                            reads = 0;
                        }
                        double dil_reads = Math.floor(reads/100.0);
                        readStats.addValue(dil_reads);

                    });

                    double mean_100 = Math.round(readStats.getMean() * 1e4) / 1e4;
                    double var_100 = Math.round(readStats.getStandardDeviation() * 1e4) / 1e4;
                    double skew_100 = Math.round(readStats.getSkewness() * 1e4) / 1e4;
                    double pooled_reads_100 = Math.round(readStats.getSum() * 1e4) / 1e4;
                    double pool_Ct_100 = ((Math.log(readStats.getSum())/Math.log(2))-37)*(-1.0);



//                    System.out.println(t_curr + "\t" + pooled_reads / (double) pool_size + "\t" + skew + "\t" + new_I / (double) N + "\t" + (Y_curr) / (double) N
//                            + "\t" + pool[0] + "\t" + pool[1] + "\t" + pool[2]);

                    try {
                        writer1.write(t_curr + "," +
                                pooled_reads_25 + "," + skew_25 + "," +
                                pooled_reads_50 + "," + skew_50 + "," +
                                pooled_reads_100 + "," + skew_100 + "," +
                                new_I / (double) N + "," + (Y_curr) / (double) N + "," +
                                X_curr + "," + Y_curr + "," + Z_curr + "," +
                                pool_25[0] + "," + pool_25[1] + "," + pool_25[2] + "," +
                                pool_50[0] + "," + pool_50[1] + "," + pool_50[2] + "," +
                                pool_100[0] + "," + pool_100[1] + "," + pool_100[2] + "," +
                                pool_Ct_25 + "," + pool_Ct_50 + "," + pool_Ct_100 + "," +
                                (ns + 1) + "\n");
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }


            }


        }
        try {
            writer1.close();
        } catch (IOException e) {
            e.printStackTrace();
        }



    }

    public void runSeasonalSIR() {


        double[] r = new double[]{0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0};
        double b = 0.09;
        double b_seas = seasonalRate(b, 1.0, 6.0, 1.0);
        double d = 0.08;
        double d_seas = seasonalRate(d, 1.0, -4.0, 2.0);
        double K = 30000.0;
        double Nr = 1000;

        int X_curr = 1000;
        int Y_curr = 0;
        int Z_curr = 0;

        double p_B = 1 - Math.exp(-1.0 * b_seas * (1 - Nr / K) * tau);
        double p_D = 1 - Math.exp(-1.0 * d_seas * tau);
        double p_SI = 1 - Math.exp(-1.0 * beta * (Y_curr / (double) Nr) * tau);
        double p_IR = 1 - Math.exp(-1.0 * gamma * tau);

        RandomEngine engine = new MersenneTwister();
        Binomial binomial = new Binomial(10, 0.2, engine);

        initialiseHistory(Y_curr, infectionHistory, currentInfected);

        int current_id = Y_curr + 1;

        double t_curr = 0.0;
        t_curr = t_curr + tau;


        while (t_curr < maxTime) {

            System.out.println(t_curr + "\t" + X_curr + "\t" + Y_curr + "\t" + Z_curr);

            int N_curr = X_curr + Y_curr + Z_curr;

            if (N_curr == 0 || p_B == 0) {
                break;
            }
            int new_X = binomial.nextInt(N_curr, p_B);
            int new_Y = 0;
            int new_Z = 0;
            int d_X = 0;
            int d_Y = 0;


            if (Y_curr > 0) {
                new_Y = binomial.nextInt(X_curr, p_SI);
                new_Z = binomial.nextInt(Y_curr, p_IR);
                d_Y = binomial.nextInt(Y_curr, p_D);

            }

            if (X_curr > 0) {
                d_X = binomial.nextInt(X_curr, p_D);
            }
            int d_Z = 0;

            if (Z_curr > 0) {
                d_Z = binomial.nextInt(Z_curr, p_D);
            }


//            if(d_Y > Y_curr) {
//                d_Y = Y_curr;
//            }
//            if(new_Z > Y_curr) {
//                new_Z = Y_curr;
//            }

            X_curr = X_curr + new_X - new_Y - d_X;
            Y_curr = Y_curr + new_Y - new_Z - d_Y;
            Z_curr = Z_curr + new_Z - d_Z;

            if (Y_curr < 0) {
                Y_curr = 0;
            }

            t_curr += tau;

            // update probabilities

            if (t_curr == 50.0) {

                Y_curr++;
                X_curr--;
                N_curr = X_curr + Y_curr + Z_curr;

            }

            b_seas = seasonalRate(b, t_curr, 6.0, 1.0);
            d_seas = seasonalRate(d, t_curr, 1.0, 1.0);


            double a = -1.0 * b_seas * (1 - N_curr / K) * tau;
            p_B = 1 - Math.exp(a);
            p_D = 1 - Math.exp(-1.0 * d_seas * tau);
            p_SI = 1 - Math.exp(-1.0 * beta * (Y_curr / (double) (N_curr)) * tau);
            p_IR = 1 - Math.exp(-1.0 * gamma * tau);


//            // new infected
//            for (int i = 0; i < new_Y; i++) {
//
//                int new_id = current_id++;
//                infectionHistory.logInfection(new_id, t_curr, -100);
//                currentInfected.add(new_id);
//            }
//
//
//            // new recovereds
//            List<Double> weights = new ArrayList();
//
//            for (Integer s : currentInfected) {
//
//                int index = infectionHistory.id.indexOf(s);
//                double time_since_infection = t_curr - infectionHistory.getBirth(index);
//                weights.add(time_since_infection);
//            }
//
//
//            List<Integer> recoveredInfections = utils.sampling.multinomSampDouble(weights, currentInfected, new_Z);
//
//            for (Integer r : recoveredInfections) {
//
//                int index = infectionHistory.id.indexOf(r);
//                infectionHistory.setDeath(index, t_curr);
//                currentInfected.remove(r);
//                recovered.add(r);
//
//            }
//
//            if (currentInfected.size() != Y_curr) {
//                System.out.println("not right " + Y_curr + "\t" + currentInfected.size());
//            }
//
//
//            // update probabilities
//
//            p_SI = 1 - Math.exp(-beta * (Y_curr / (double) N));
//            p_IR = 1 - Math.exp(-gamma);
//
//
//            t_curr += tau;
//
//            if (Y_curr == 0) {
//                break;
//            }
//
//            if (t_curr % 5 == 0) {
//
//                int[] pool = sampleIndividuals(X_curr, Y_curr, Z_curr);
//
//                List<Double> CtDist = getCtDistribution(pool, t_curr);
//
//                ctStats.clear();
//                readStats.clear();
//                CtDist.forEach(Ct -> {
////                    if(Ct < 37.0) {
////                        Ct = Ct + poolDilutionFactor;
////                    }
////                    else{
////                        Ct = 40.0;
////                    }
//                    ctStats.addValue(Ct);
//
//                    double reads = Math.round(Math.pow(2, -0.47 * Ct + 14.4));
//                    if (reads < 0) {
//                        reads = 0;
//                    }
//                    readStats.addValue(reads);
//                });
//
//                double mean = Math.round(ctStats.getMean() * 1e4) / 1e4;
//                double var = Math.round(ctStats.getStandardDeviation() * 1e4) / 1e4;
//                double skew = Math.round(ctStats.getSkewness() * 1e4) / 1e4;
//                double pooled_reads = Math.round(readStats.getSum() * 1e4) / 1e4;
//
//                System.out.println(t_curr + "\t" + pooled_reads / (double) pool_size + "\t" + skew + "\t" + new_I / (double) N + "\t" + (Y_curr) / (double) N
//                        + "\t" + pool[0] + "\t" + pool[1] + "\t" + pool[2]);


        }


    }

    public int[] sampleIndividuals(int X_curr, int Y_curr, int Z_curr, int pool_size) {

        double p_x = (double) X_curr / (double) N;
        double p_y = (double) Y_curr / (double) N;
        double p_z = (double) Z_curr / (double) N;

        int sample_X = Binomial.staticNextInt(pool_size, p_x);
        int sample_Y = Binomial.staticNextInt(pool_size, p_y);

        if ((sample_X + sample_Y) > pool_size) {
            sample_Y = 0;
        }
        int sample_Z = pool_size - sample_X - sample_Y;

        return new int[]{sample_X, sample_Y, sample_Z};
    }

    public int[] sampleIndividuals(int X_curr, int E_curr, int Y_curr, int Z_curr, int pool_size) {

        double p_x = (double) X_curr / (double) N;
        double p_ey = (double) (E_curr + Y_curr) / (double) N;
        double p_z = (double) Z_curr / (double) N;

        int sample_X = Binomial.staticNextInt(pool_size, p_x);
        int sample_EY = Binomial.staticNextInt(pool_size, p_ey);

        if ((sample_X + sample_EY) > pool_size) {
            sample_EY = 0;
        }
        int sample_Z = pool_size - sample_X - sample_EY;

        return new int[]{sample_X, sample_EY, sample_Z};
    }


    public List<Double> getCtDistribution(int[] pool, double t_curr) {

        List<Integer> exposed_and_infected = new ArrayList<>();
        exposed_and_infected.addAll(exposed);
        exposed_and_infected.addAll(currentInfected);

        Double[] Ct_susceptible = new Double[pool[0]];
        Arrays.fill(Ct_susceptible, 40.0);
        Double[] Ct_infected = getSampleCt(pool[1], exposed_and_infected, t_curr);
        Double[] Ct_recovered = getSampleCt(pool[2], recovered, t_curr);


        List<Double> Ct_all = new ArrayList<>();

        Ct_all.addAll(Arrays.asList(Ct_susceptible));
        Ct_all.addAll(Arrays.asList(Ct_infected));
        Ct_all.addAll(Arrays.asList(Ct_recovered));

        return Ct_all;
    }


    public Double[] getSampleCt(int sample_size, List<Integer> infections, double t_curr) {

        if (sample_size == 0) {
            Double[] doubles = new Double[0];
            return doubles;
        } else {
            Double[] Ct = new Double[sample_size];

            List<Integer> currentInfected_temp = new ArrayList<>();
            currentInfected_temp.addAll(infections);
            Collections.shuffle(currentInfected_temp);

            for (int i = 0; i < sample_size; i++) {

                Integer infection = currentInfected_temp.get(i);
                int index = infectionHistory.id.indexOf(infection);
                double t_birth = infectionHistory.getBirth(index);

                double t = t_curr - t_birth;
                double Ct_value = CtModel(t);

                Ct[i] = Ct_value;


            }

            return Ct;
        }


    }

    public double seasonalRate(double b, double t, double freq, double amp) {

        double rate = b * (amp + Math.sin(2 * Math.PI * tau * t / freq));


        return rate;

    }

    public double CtModel(double t) {

        // see Hay et al Science (2021) Supp material Ct Value model
        double t_eclipse = 0; // (0 days) Time from infection to initial viral growth
        double t_peak = 4; // (5 days ) Time from initial viral growth to peak viral load
        double t_switch = 8; // (9.38 days) Time from peak viral load to secondary waning phase
        double t_mod = 14; // (14 days) time from secondary waning phase until gumbel distribution reaches its min scale parameter
        double t_LOD = Double.POSITIVE_INFINITY; // ( inf days ) Time from infection until modal Ct value is equal to the limit of detection

        double sigma_obs = 5; // Initial scale parameter for the Gumbel distribution until a=teclipse+tpeak+tswitch
        double s_mod = 0.4; // 0.4 multiplicative factor applied to scale paramter for the Gumble distrbution - starting at t_eclipse + t_peak + t_switch + t_scle
        double C_zero = 40; // Ct value at time of infection
        double C_peak = 20; // (20) Modal Ct value at peak viral load
        double C_switch = 30; // (33) Modal Ct value at a = teclipse + tpeak + tswitch
        double C_LOD = 36; // Limit of detection of Ct value

        double C_mode_t = -1;
        double sigma_t = -1;


        if (t <= t_eclipse) {

            C_mode_t = C_zero;
        } else if ((t_eclipse < t) && (t <= (t_eclipse + t_peak))) {

            C_mode_t = C_zero + ((C_peak - C_zero) / (t_peak)) * (t - t_eclipse);

        } else if (((t_eclipse + t_peak) < t) && (t <= (t_eclipse + t_peak + t_switch))) {

            C_mode_t = C_peak + ((C_switch - C_peak) / t_switch) * (t - t_eclipse - t_peak);
        } else if ((t_eclipse + t_peak + t_switch) < t) {

            C_mode_t = C_switch + ((C_LOD - C_switch) / (t_LOD - t_switch - t_peak - t_eclipse)) * (t - t_eclipse - t_peak - t_switch);

        }

        if (t < (t_eclipse + t_peak + t_switch)) {

            sigma_t = sigma_obs;

        } else if (((t_eclipse + t_peak + t_switch) <= t) && (t < (t_eclipse + t_peak + t_switch + t_mod))) {

            sigma_t = sigma_obs * (1 - ((1 - s_mod) / t_mod) * (t - t_eclipse - t_peak - t_switch));

        } else if (((t_eclipse + t_peak + t_mod) <= t)) {

            sigma_t = sigma_obs * s_mod;

        }

        // drawing from a Gumbel dist Ct ~ (C_mode_t, sigma_t)
        double Ct = C_mode_t - sigma_t * Math.log(Exponential.staticNextDouble(1));

        return Ct;

    }

    public void initialiseHistory(int Y_curr, infections history, List currentInfected) {

        for (int i = 0; i < Y_curr; i++) {
            history.logInfection(i + 1, 0, -100);
            currentInfected.add(i + 1);

        }
    }

    public static void main(String[] args) {

        stochasticSIR test = new stochasticSIR();

        test.runSIR(100);

//        DescriptiveStatistics stats = new DescriptiveStatistics();
//        for(int i = 0; i < 100; i ++) {
//          double Ct =  test.CtModel(14);
//          stats.addValue(Ct);
//        }
//        System.out.println(stats.getMean());
//        System.out.println(Math.pow(stats.getVariance(), 0.5));
//        System.out.println(stats.getPercentile(25) + "\t" + stats.getPercentile(50) + "\t" +stats.getPercentile(75));


//        test.runSIR();

    }
}
