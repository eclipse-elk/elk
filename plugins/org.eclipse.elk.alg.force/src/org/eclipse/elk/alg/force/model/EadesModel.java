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
import org.eclipse.elk.alg.force.graph.FParticle;
import org.eclipse.elk.alg.force.options.ForceOptions;
import org.eclipse.elk.core.math.KVector;

/**
 * A force model after the Eades approach.
 *
 * @author msp
 */
public final class EadesModel extends AbstractForceModel {

    /** factor used for force calculations when the distance of two particles is zero. */
    private static final double ZERO_FACTOR = 100;

    /** the maximal number of iterations after which the model stops. */
    private int maxIterations = ForceOptions.ITERATIONS.getDefault();
    /** the spring length that determines the optimal distance of connected nodes. */
    private double springLength = ForceOptions.SPACING_NODE_NODE.getDefault();
    /** additional factor for repulsive forces. */
    private double repulsionFactor = ForceOptions.REPULSION.getDefault();
    
    @Override
    protected void initialize(final FGraph graph) {
        super.initialize(graph);
        maxIterations = graph.getProperty(ForceOptions.ITERATIONS);
        springLength = graph.getProperty(ForceOptions.SPACING_NODE_NODE).doubleValue();
        repulsionFactor = graph.getProperty(ForceOptions.REPULSION).doubleValue();
    }
    
    @Override
    protected boolean moreIterations(final int count) {
        return count < maxIterations;
    }

    @Override
    protected KVector calcDisplacement(final FParticle forcer, final FParticle forcee) {
        avoidSamePosition(getRandom(), forcer, forcee);

        // compute distance (z in the original algorithm)
        KVector displacement = forcee.getPosition().clone().sub(forcer.getPosition());
        double length = displacement.length();
        double d = Math.max(0, length - forcer.getRadius() - forcee.getRadius());
        
        // calculate attractive or repulsive force, depending of adjacency
        double force;
        int connection = getGraph().getConnection(forcer, forcee);
        if (connection > 0) {
            force = -attractive(d, springLength) * connection;
        } else {
            force = repulsive(d, repulsionFactor) * forcer.getProperty(ForceOptions.PRIORITY);
        }

        // scale distance vector to the amount of repulsive forces
        displacement.scale(force / length);

        return displacement;
    }
    
    /**
     * Compute repulsion force between the forcee and the forcer.
     *
     * @param d the distance between the two particles
     * @param r the factor for repulsive force
     * @return a force exerted on the forcee 
     */
    private static double repulsive(final double d, final double r) {
        if (d > 0) {
            return r / (d * d);
        } else {
            return r * ZERO_FACTOR;
        }
    }
    
    /**
     * Compute attraction force between the forcee and the forcer.
     * 
     * @param d the distance between the two particles
     * @param s the spring length
     * @return a force exerted on the forcee
     */
    public static double attractive(final double d, final double s) {
        if (d > 0) {
            return Math.log(d / s);
        } else {
            return -ZERO_FACTOR;
        }
    }

}
