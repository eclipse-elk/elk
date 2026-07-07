/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.force.model;

import java.util.List;
import java.util.Random;

import org.eclipse.elk.alg.force.graph.FBendpoint;
import org.eclipse.elk.alg.force.graph.FEdge;
import org.eclipse.elk.alg.force.graph.FGraph;
import org.eclipse.elk.alg.force.graph.FLabel;
import org.eclipse.elk.alg.force.graph.FNode;
import org.eclipse.elk.alg.force.graph.FParticle;
import org.eclipse.elk.alg.force.options.ForceOptions;
import org.eclipse.elk.alg.force.options.InternalProperties;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Superclass for force models.
 *
 * @author msp
 */
public abstract class AbstractForceModel {

    /** the random number generator used for this model. */
    private Random random;
    /** the graph for which this model was initialized. */
    private FGraph graph;
    /** upper bound for particle coordinates. */
    private double dispBound;
    
    /** factor by which nodes influence the displacement bound. */
    private static final double DISP_BOUND_FACTOR = 16;
    
    /**
     * Initialize the force model with the given graph. Subclasses that override this
     * must call the superclass method first.
     * 
     * @param fgraph a force graph
     */
    protected void initialize(final FGraph fgraph) {
        this.graph = fgraph;
        this.random = fgraph.getProperty(InternalProperties.RANDOM);
        
        // calculate the adjacency matrix for the graph
        fgraph.calcAdjacency();
        
        // calculate an upper bound for particle displacement
        dispBound = Math.max(fgraph.getNodes().size() * DISP_BOUND_FACTOR + fgraph.getEdges().size(),
                DISP_BOUND_FACTOR * DISP_BOUND_FACTOR);
        
        // if interactive mode is off, randomize the layout
        if (!fgraph.getProperty(ForceOptions.INTERACTIVE)) {
            double posScale = graph.getNodes().size();
            for (FNode node : fgraph.getNodes()) {
                KVector pos = node.getPosition();
                pos.x = random.nextDouble() * posScale;
                pos.y = random.nextDouble() * posScale;
            }
        }
        
        // create bend points for node repulsion
        List<FBendpoint> bends = fgraph.getBendpoints();
        for (FEdge edge : fgraph.getEdges()) {
            int count = edge.getProperty(ForceOptions.REPULSIVE_POWER);
            if (count > 0) {
                for (int i = 0; i < count; i++) {
                    bends.add(new FBendpoint(edge));
                }
                edge.distributeBendpoints();
            }
        }
    }
    
    /**
     * Return the graph for which the model was initialized.
     * 
     * @return the force graph
     */
    public final FGraph getGraph() {
        return graph;
    }
    
    /**
     * Return the random number generator for this model.
     * 
     * @return the random number generator
     */
    protected final Random getRandom() {
        return random;
    }
    
    /**
     * Perform layout on the given force graph.
     * 
     * @param fgraph a force graph
     * @param monitor a progress monitor
     */
    public void layout(final FGraph fgraph, final IElkProgressMonitor monitor) {
        monitor.begin("Component Layout", 1);
        initialize(fgraph);
        int iterations = 0;
        
        while (moreIterations(iterations) && !monitor.isCanceled()) {

            iterationDone();
            // calculate attractive and repulsive forces
            for (FParticle v : fgraph.getParticles()) {
                for (FParticle u : fgraph.getParticles()) {
                    if (u != v) {
                        KVector displacement = calcDisplacement(u, v);
                        if (displacement != null) {
                            v.getDisplacement().add(displacement);
                        }
                    }
                }
            }
            
            // apply calculated displacement
            for (FParticle v : fgraph.getParticles()) {
                KVector d = v.getDisplacement();
                d.bound(-dispBound, -dispBound, dispBound, dispBound);
                v.getPosition().add(d);
                d.reset();
            }
            iterations++;
        }
        monitor.done();
    }
    
    /**
     * Perform all necessary calculations after a full iteration. Subclasses must call
     * the superclass method first.
     */
    protected void iterationDone() {
        for (FEdge edge : graph.getEdges()) {
            // adjust label positions
            for (FLabel label : edge.getLabels()) {
                label.refreshPosition();
            }
            
            // adjust bend point positions
            edge.distributeBendpoints();
        }
    }
    
    /**
     * Determine whether more iterations need to be performed.
     * 
     * @param count the current iterations count
     * @return {@code true} if more iterations shall be done
     */
    protected abstract boolean moreIterations(int count);
    
    /**
     * Calculate a displacement for the given particles.
     * 
     * @param forcer the particle that is causing the force
     * @param forcee the particle that is affected by the force
     * @return a displacement vector for the forcee, or {@code null} if no force is applied
     */
    protected abstract KVector calcDisplacement(FParticle forcer, FParticle forcee);
    
    /**
     * Avoid having nodes on the same position by moving them a little.
     * 
     * @param random the random number generator
     * @param u first node
     * @param v second node
     */
    public static void avoidSamePosition(final Random random, final FParticle u, final FParticle v) {
        KVector pu = u.getPosition();
        KVector pv = v.getPosition();
        while (pu.x - pv.x == 0 && pu.y - pv.y == 0) {
            boolean triedForBendPoints = false;
            if (u instanceof FBendpoint && v instanceof FBendpoint && !triedForBendPoints) {
                // Wiggle orthogonal to edge direction
                FEdge uE = ((FBendpoint) u).getEdge();
                KVector uVector = new KVector(uE.getTargetPoint()).sub(uE.getSourcePoint());
                double length = 2;
                KVector orthogonaluV = new KVector(uVector.x / uVector.length() * length,
                        -uVector.y / uVector.length() * length);
                pu.add(orthogonaluV);

                FEdge vE = ((FBendpoint) v).getEdge();
                KVector vVector = new KVector(vE.getTargetPoint()).sub(vE.getSourcePoint());
                length = uVector == vVector ? -2 : 2;
                KVector orthogonalvV = new KVector((vVector.x / vVector.length()) * length,
                        -(vVector.y / vVector.length()) * length);
                pu.add(orthogonalvV);
                triedForBendPoints = true;
                
            } else {
                pu.wiggle(random, 1);
                pv.wiggle(random, 1);
            }
        }
    }

}
