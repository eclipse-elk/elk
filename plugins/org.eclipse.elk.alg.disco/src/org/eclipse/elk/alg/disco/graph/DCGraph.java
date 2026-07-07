/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.disco.graph;

import java.util.Collection;
import java.util.Set;

import org.eclipse.elk.alg.disco.ICompactor;
import org.eclipse.elk.alg.disco.transform.IGraphTransformer;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.graph.properties.MapPropertyHolder;

import com.google.common.collect.Sets;

/**
 * This class represents the underlying graph data structure for the compaction of connected components. The benefit of
 * using an intermediate graph structure is that all classes implementing {@link ICompactor} can be used with all
 * possible graph data structures, as long as an {@link IGraphTransformer} is implemented for them.
 */
public class DCGraph extends MapPropertyHolder {
    /** the serial version UID. */
    private static final long serialVersionUID = 2230448069992361373L;

    ///////////////////////////////////////////////////////////////////////////////
    // Variables

    /** All the different connected {@link DCComponent DCComponents} this graph consists of. */
    private Set<DCComponent> components = Sets.newLinkedHashSet();
    /** Width and height of this graph (after compaction). */
    private KVector dimensions;
    /** analogous to the insets from the (old) KGraph. */
    private double inset = 0.0f;

    /**
     * Takes collections of {@link DCElement DCElements} nested in an outer collection, interprets each inner collection
     * as a connected component of the {@link DCGraph} and constructs the graph accordingly; including instantiation of
     * each DCComponent and setting the parent component of each {@link DCElement}.
     * 
     * @param components
     *            Collection of future {@link DCComponent DCComponents} given as collections of {@link DCElement
     *            DCElements}
     */
    public <C extends Collection<DCElement>, D extends Collection<C>> DCGraph(final D components) {
        for (C elements : components) {
            DCComponent component = new DCComponent();
            component.addElements(elements);
            this.components.add(component);
        }
    }

    /**
     * Takes collections of {@link DCElement DCElements} nested in an outer collection, interprets each inner collection
     * as a connected component of the {@link DCGraph} and constructs the graph accordingly; including instantiation of
     * each DCComponent and setting the parent component of each {@link DCElement}.
     * 
     * @param components
     *            Collection of future {@link DCComponent DCComponents} given as collections of {@link DCElement
     *            DCElements}
     * @param inset
     *            optional insets
     */
    public <C extends Collection<DCElement>, D extends Collection<C>> DCGraph(final D components, final double inset) {
        this(components);
        this.inset = inset;
    }

    /**
     * Gets all {@link DCComponent DCComponents} of this graph.
     * 
     * @return {@link DCComponent DCComponents} of this graph
     */
    public Set<DCComponent> getComponents() {
        return components;
    }

    /**
     * Returns width and height of the graph structure.
     * 
     * @return the dimensions
     */
    public KVector getDimensions() {
        return dimensions;
    }

    /**
     * Sets the width and height of the graph structure (called after compaction).
     * 
     * @param dimensions
     *            the dimensions to set
     */
    public void setDimensions(final KVector dimensions) {
        this.dimensions = dimensions;
    }

    /**
     * Returns the inset.
     * 
     * @return inset
     */
    public double getInset() {
        return inset;
    }

}
