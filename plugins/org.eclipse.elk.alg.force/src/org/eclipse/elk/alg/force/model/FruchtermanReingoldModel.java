/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.force.model;

import org.eclipse.elk.alg.force.graph.FGraph;
import org.eclipse.elk.alg.force.graph.FNode;
import org.eclipse.elk.alg.force.graph.FParticle;
import org.eclipse.elk.alg.force.options.ForceOptions;
import org.eclipse.elk.core.math.KVector;

/**
 * A force model after the Fruchterman-Reingold approach.
 *
 * @author msp
 * @author tmn
 * @author fhol
 */
public final class FruchtermanReingoldModel extends AbstractForceModel {

    /** factor that determines the C constant used for calculation of K. */
    private static final double SPACING_FACTOR = 0.01;
    /** factor used for repulsive force when the distance of two particles is zero. */
    private static final double ZERO_FACTOR = 100;
    
    /** the current temperature of the system. */
    private double temperature = ForceOptions.TEMPERATURE.getDefault();
    /** the temperature threshold for stopping the model. */
    private double threshold;
    /** the main constant used for force calculations. */
    private double k;
    
    @Override
    protected void initialize(final FGraph graph) {
        super.initialize(graph);
        temperature = graph.getProperty(ForceOptions.TEMPERATURE).doubleValue();
        threshold = temperature / graph.getProperty(ForceOptions.ITERATIONS);
        
        // calculate an appropriate value for K
        int n = graph.getNodes().size();
        double totalWidth = 0;
        double totalHeight = 0;
        for (FNode v : graph.getNodes()) {
            totalWidth += v.getSize().x;
            totalHeight += v.getSize().y;
        }
        double area = totalWidth * totalHeight;
        double c = graph.getProperty(ForceOptions.SPACING_NODE_NODE) * SPACING_FACTOR;
        k = Math.sqrt(area / (2 * n)) * c;
    }

    @Override
    protected boolean moreIterations(final int count) {
        return temperature > 0;
    }

    @Override
    protected KVector calcDisplacement(final FParticle forcer, final FParticle forcee) {
        avoidSamePosition(getRandom(), forcer, forcee);

        // compute distance (z in the original algorithm)
        KVector displacement = forcee.getPosition().clone().sub(forcer.getPosition());
        double length = displacement.length();
        double d = Math.max(0, length - forcer.getRadius() - forcee.getRadius());
        
        // calculate repulsive force, independent of adjacency
        double force = repulsive(d, k) * forcer.getProperty(ForceOptions.PRIORITY);
        
        // calculate attractive force, depending of adjacency
        int connection = getGraph().getConnection(forcer, forcee);
        if (connection > 0) {
            force -= attractive(d, k) * connection;
        }

        // scale distance vector to the amount of repulsive forces
        displacement.scale(force * temperature / length);

        return displacement;
    }
    
    @Override
    protected void iterationDone() {
        super.iterationDone();
        temperature -= threshold;
    }
    
    /**
     * Calculates the amount of repulsive force along a distance.
     * 
     * @param d the distance over which the force exerts
     * @param k the space parameter, depending on the available area
     * @return the amount of the repulsive force
     */
    private static double repulsive(final double d, final double k) {
        if (d > 0) {
            return k * k / d;
        } else {
            return k * k * ZERO_FACTOR;
        }
    }
    
    /**
     * Calculates the amount of attracting force along a distance.
     * 
     * @param d the distance over which the force exerts
     * @param k the space-parameter, depending on the available area
     * @return the amount of the attracting force
     */
    public static double attractive(final double d, final double k) {
        return d * d / k;
    }

}
