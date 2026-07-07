/*******************************************************************************
 * Copyright (c) 2010, 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.graph;

import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

import org.eclipse.elk.core.math.KVector;

import com.google.common.collect.Lists;

/**
 * A layered graph has a set of layers that contain the nodes, as well as a
 * list of nodes that are not yet assigned to a layer. Layout algorithms are
 * required to layout the graph from left to right. If another layout direction
 * is desired, it can be obtained by pre-processing and post-processing the graph.
 * In contrast to the KGraph structure, the LGraph is not EMF-based, but plain Java.
 * It is optimized for being processed by a layer-based layout algorithm.
 */
public final class LGraph extends LGraphElement implements Iterable<Layer> {
    
    /**
     * The serial version UID.
     */
    private static final long serialVersionUID = -8006835373897072852L;
    
    /**
     * The total size of the drawing. The total size does not include padding. If the
     * actual size, including padding, is required, this can be obtained by calling
     * {@link #getActualSize()}.
     */
    private final KVector size = new KVector();
    
    /**
     * The graph's effective padding. This already includes the padding set on the imported graph by the
     * client, as well as any space required for inside node labels. These are the effective paddings
     * returned to the client if {@link org.eclipse.elk.core.options.SizeOptions#COMPUTE_PADDING
     * SizeOptions.COMPUTE_PADDING}
     * is active.
     */
    private final LPadding padding = new LPadding();
    
    /**
     * The offset that will be applied to all node positions. Adding the offset to the coordinates of
     * child nodes results in child node coordinates that are relative to the node's client area (which
     * does not include the padding around it).
     */
    private final KVector offset = new KVector();
    
    /**
     * Nodes that are not currently part of a layer.
     */
    private final List<LNode> layerlessNodes = Lists.newArrayList();
    
    /**
     * The layers of the layered graph.
     */
    private final List<Layer> layers = Lists.newArrayList();
    
    /**
     * The parent node in which this graph is nested, or {@code null}.
     */
    private LNode parentNode;
    
    
    /**
     * Returns the size of the graph, that is the bounding box that covers the
     * whole drawing. The size does not include padding. Modifying the
     * returned value changes the size of the graph.
     * 
     * @return the size of the layered graph; modify to change the graph size.
     */
    public KVector getSize() {
        return size;
    }
    
    /**
     * Returns the graph's size including any borders. If the graph represents a
     * hierarchical node, the returned size represents the node's size. The returned
     * size can be modified at will without having any influence on the graph's size
     * or the actual size returned on the next method call.
     * 
     * @return the graph's size including borders.
     */
    public KVector getActualSize() {
        return new KVector(
                size.x + padding.left + padding.right,
                size.y + padding.top + padding.bottom);
    }
    
    /**
     * Returns the padding of the graph. The padding determines the amount of space between
     * the content area and the graph's actual border. Modifying the returned value
     * changes the padding.
     * 
     * @return the padding; modify to change the graph's padding.
     */
    public LPadding getPadding() {
        return padding;
    }

    /**
     * Returns the offset for the graph, that is a coordinate vector that has
     * to be added to all position values of nodes and edges. It is usually used
     * to reserve some space in the content area for additional edge routing.
     * 
     * <b>Note:</b> Since many different parts of the algorithm may contribute to
     * the offset, never set the offset to an absolute value! Rather, only add to
     * the offset!
     * 
     * @return the offset of the layered graph
     */
    public KVector getOffset() {
        return offset;
    }
    
    /**
     * Returns the list of nodes that are not currently assigned to a layer.
     * When creating a graph, put the nodes here.
     * 
     * @return the layerless nodes.
     */
    public List<LNode> getLayerlessNodes() {
        return layerlessNodes;
    }

    /**
     * Returns the list of layers of the graph. Layers are created automatically by the layout
     * algorithm, so this list must not be touched when the graph is created.
     * 
     * @return the layers
     */
    public List<Layer> getLayers() {
        return layers;
    }
    
    /**
     * Returns the parent node in which this graph is nested, if any.
     * 
     * @return the parent node or {@code null}
     */
    public LNode getParentNode() {
        return parentNode;
    }
    
    /**
     * Sets the parent node in which this graph is nested.
     * 
     * @param parentNode the new parent node
     */
    public void setParentNode(final LNode parentNode) {
        this.parentNode = parentNode;
    }

    /**
     * Returns an iterator over the layers.
     * 
     * @return an iterator for the layers of this layered graph
     */
    public Iterator<Layer> iterator() {
        return layers.iterator();
    }
    
    /**
     * Returns a two-dimensional array representation of all nodes in the graph. The first dimension
     * represents the layers, the second dimension represents the nodes in each layer.
     * 
     * @return two-dimensional array of nodes in the graph.
     */
    public LNode[][] toNodeArray() {
        LNode[][] lgraphArray = new LNode[layers.size()][];
        ListIterator<Layer> layerIter = layers.listIterator();
        while (layerIter.hasNext()) {
            Layer layer = layerIter.next();
            int layerIndex = layerIter.previousIndex();
            lgraphArray[layerIndex] = LGraphUtil.toNodeArray(layer.getNodes());
        }
        
        return lgraphArray;
    }
    
    @Override
    public String toString() {
        if (layers.isEmpty()) {
            return "G-unlayered" + layerlessNodes.toString();
        } else if (layerlessNodes.isEmpty()) {
            return "G-layered" + layers.toString();
        }
        return "G[layerless" + layerlessNodes.toString() + ", layers" + layers.toString() + "]";
    }
    
}
