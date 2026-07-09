/*******************************************************************************
 * Copyright (c) 2017, 2020 Kiel University and others.
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
package org.eclipse.elk.alg.common.nodespacing.internal;

import java.util.EnumMap;
import java.util.Set;

import org.eclipse.elk.alg.common.nodespacing.cellsystem.AtomicCell;
import org.eclipse.elk.alg.common.nodespacing.cellsystem.ContainerArea;
import org.eclipse.elk.alg.common.nodespacing.cellsystem.GridContainerCell;
import org.eclipse.elk.alg.common.nodespacing.cellsystem.LabelCell;
import org.eclipse.elk.alg.common.nodespacing.cellsystem.StripContainerCell;
import org.eclipse.elk.alg.common.nodespacing.cellsystem.StripContainerCell.Strip;
import org.eclipse.elk.core.UnsupportedConfigurationException;
import org.eclipse.elk.core.math.ElkMargin;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.NodeLabelPlacement;
import org.eclipse.elk.core.options.PortAlignment;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortLabelPlacement;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.options.SizeConstraint;
import org.eclipse.elk.core.options.SizeOptions;
import org.eclipse.elk.core.util.IndividualSpacings;
import org.eclipse.elk.core.util.adapters.GraphAdapters.GraphAdapter;
import org.eclipse.elk.core.util.adapters.GraphAdapters.NodeAdapter;

import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.TreeMultimap;

/**
 * Data holder class to be passed around to avoid having too much state in the size calculation classes. Some of the
 * most relevant settings are copied into variables for convenience. The node's size is stored in a separate variable
 * in this context, to be used by the algorithm. Once the algorithm has finished it can apply the calculated size to
 * the node by calling {@link #applyNodeSize()}.
 */
public final class NodeContext {
   
    // CHECKSTYLEOFF Visibility
    // This is a purely internal data holder class, so we ditch setters for public fields.
    
    /////////////////////////////////////////////////////////////////////////////////
    // Convenience Access to Things
    
    /** The node we calculate stuff for. */
    public final NodeAdapter<?> node;
    /** The node's size. This will be used during the algorithm, to be applied (or not) once it is finished. */
    public final KVector nodeSize;
    /** Whether this node has stuff inside it or not. */
    public final boolean treatAsCompoundNode;
    /** The node's size constraints. */
    public final Set<SizeConstraint> sizeConstraints;
    /** The node's size options. */
    public final Set<SizeOptions> sizeOptions;
    /** Port constraints set on the node. */
    public final PortConstraints portConstraints;
    /** Whether port labels are placed inside or outside. */
    public final Set<PortLabelPlacement> portLabelsPlacement;
    /** Whether to treat port labels as a group when centering them next to eastern or western ports. */
    public final boolean portLabelsTreatAsGroup;
    /** Where node labels are placed by default. */
    public final Set<NodeLabelPlacement> nodeLabelPlacement;
    /** Space to leave around the node label area. */
    public final ElkPadding nodeLabelsPadding;
    /** Space between a node and its outside labels. */
    public final double nodeLabelSpacing;
    /** Space between two labels. */
    public final double labelLabelSpacing;
    /** Space between two different label cells. */
    public final double labelCellSpacing;
    /** Space between a port and another port. */
    public final double portPortSpacing;
    /** Horizontal space between a port and its labels. */
    public final double portLabelSpacingHorizontal;
    /** Vertical space between a port and its labels. */
    public final double portLabelSpacingVertical;
    /** Margin to leave around the set of ports on each side. */
    public final ElkMargin surroundingPortMargins;
    /** Whether node is being laid out in top-down layout mode. */
    public final boolean topdownLayout;


    /////////////////////////////////////////////////////////////////////////////////
    // More Contexts
    
    /** Context objects that hold more information about each port. Sorted left-to-right / top-to-bottom. */
    public final Multimap<PortSide, PortContext> portContexts = TreeMultimap.create(
            NodeContext::comparePortSides, NodeContext::comparePortContexts);
    
    
    /////////////////////////////////////////////////////////////////////////////////
    // The Cell System
    
    /** The main cell that holds all the cells that make up the node. */
    public final StripContainerCell nodeContainer;
    /** The main cell's middle row, which will contain further cells. */
    public final StripContainerCell nodeContainerMiddleRow;
    /** The grid container that represents the node's area reserved for inside node labels (and the client area). */
    public GridContainerCell insideNodeLabelContainer;
    /**
     * All cells that will describe the space required for ports and for inside port labels. The paddings on these
     * things not only describe the space to be left between ports and their inside labels, but also the space to be
     * left between ports and the node border (the additional port spacing). These paddings will be updated as we
     * calculate more information over the course of the algorithm.
     */
    public final EnumMap<PortSide, AtomicCell> insidePortLabelCells = Maps.newEnumMap(PortSide.class);
    /** All container cells that will hold label cells for outside node labels. */
    public final EnumMap<PortSide, StripContainerCell> outsideNodeLabelContainers = Maps.newEnumMap(PortSide.class);
    /** All of the label cells created for possible node labels, both inside and outside. */
    public final EnumMap<NodeLabelLocation, LabelCell> nodeLabelCells = Maps.newEnumMap(NodeLabelLocation.class);
    
    
    /////////////////////////////////////////////////////////////////////////////////
    // Creation
    
    /**
     * Creates a new context object for the given node, fully initialized with the node's settings.
     * 
     * @param parentGraph the node's parent graph.
     * @param node the node to create the context for.
     */
    public NodeContext(final GraphAdapter<?> parentGraph, final NodeAdapter<?> node) {
        this.node = node;
        this.nodeSize = new KVector(node.getSize());
        
        // Top-down layout
        topdownLayout = node.getProperty(CoreOptions.TOPDOWN_LAYOUT);
        
        // Compound node
        treatAsCompoundNode = node.isCompoundNode() || node.getProperty(CoreOptions.INSIDE_SELF_LOOPS_ACTIVATE);
        
        // Core size settings
        sizeConstraints = node.getProperty(CoreOptions.NODE_SIZE_CONSTRAINTS);
        sizeOptions = node.getProperty(CoreOptions.NODE_SIZE_OPTIONS);
        portConstraints = node.getProperty(CoreOptions.PORT_CONSTRAINTS);
        portLabelsPlacement = node.getProperty(CoreOptions.PORT_LABELS_PLACEMENT);
        if (!PortLabelPlacement.isValid(portLabelsPlacement)) {
            throw new UnsupportedConfigurationException("Invalid port label placement: " + portLabelsPlacement);
        }

        portLabelsTreatAsGroup = node.getProperty(CoreOptions.PORT_LABELS_TREAT_AS_GROUP);
        nodeLabelPlacement = node.getProperty(CoreOptions.NODE_LABELS_PLACEMENT);
        if (!NodeLabelPlacement.isValid(nodeLabelPlacement)) {
            throw new UnsupportedConfigurationException("Invalid node label placement: " + nodeLabelPlacement);
        }
        
        // Copy spacings for convenience
        nodeLabelsPadding = IndividualSpacings.getIndividualOrInherited(node, CoreOptions.NODE_LABELS_PADDING);
        nodeLabelSpacing = IndividualSpacings.getIndividualOrInherited(node, CoreOptions.SPACING_LABEL_NODE);
        labelLabelSpacing = IndividualSpacings.getIndividualOrInherited(node, CoreOptions.SPACING_LABEL_LABEL);
        portPortSpacing = IndividualSpacings.getIndividualOrInherited(node, CoreOptions.SPACING_PORT_PORT);
        portLabelSpacingHorizontal =
                IndividualSpacings.getIndividualOrInherited(node, CoreOptions.SPACING_LABEL_PORT_HORIZONTAL);
        portLabelSpacingVertical =
                IndividualSpacings.getIndividualOrInherited(node, CoreOptions.SPACING_LABEL_PORT_VERTICAL);
        surroundingPortMargins = IndividualSpacings.getIndividualOrInherited(
                node, CoreOptions.SPACING_PORTS_SURROUNDING);
        
        labelCellSpacing = 2 * labelLabelSpacing;
        
        // Create main cells (the others will be created later)
        boolean symmetry = !sizeOptions.contains(SizeOptions.ASYMMETRICAL);
        nodeContainer = new StripContainerCell(Strip.VERTICAL, symmetry, 0);
        
        nodeContainerMiddleRow = new StripContainerCell(Strip.HORIZONTAL, symmetry, 0);
        nodeContainer.setCell(ContainerArea.CENTER, nodeContainerMiddleRow);
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////
    // Application
    
    /**
     * Applies the node size stored in this context to the actual node.
     */
    public void applyNodeSize() {
        node.setSize(nodeSize);
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////
    // Utility Methods
    
    /**
     * Returns the port alignment that applies to the given side of the node.
     * 
     * @param portSide
     *            the side.
     * @return the side's port alignment.
     */
    public PortAlignment getPortAlignment(final PortSide portSide) {
        PortAlignment alignment = null;
        
        switch (portSide) {
        case NORTH:
            if (node.hasProperty(CoreOptions.PORT_ALIGNMENT_NORTH)) {
                alignment = node.getProperty(CoreOptions.PORT_ALIGNMENT_NORTH);
            }
            break;

        case SOUTH:
            if (node.hasProperty(CoreOptions.PORT_ALIGNMENT_SOUTH)) {
                alignment = node.getProperty(CoreOptions.PORT_ALIGNMENT_SOUTH);
            }
            break;

        case EAST:
            if (node.hasProperty(CoreOptions.PORT_ALIGNMENT_EAST)) {
                alignment = node.getProperty(CoreOptions.PORT_ALIGNMENT_EAST);
            }
            break;

        case WEST:
            if (node.hasProperty(CoreOptions.PORT_ALIGNMENT_WEST)) {
                alignment = node.getProperty(CoreOptions.PORT_ALIGNMENT_WEST);
            }
            break;
        }
        
        // Fall back to basic port alignment if we haven't found a more specific one yet
        if (alignment == null) {
            alignment = node.getProperty(CoreOptions.PORT_ALIGNMENT_DEFAULT);
        }
        
        return alignment;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////
    // Comparators
    
    /**
     * Comparator for two port sides.
     */
    public static int comparePortSides(final PortSide portSide1, final PortSide portSide2) {
        return Integer.compare(portSide1.ordinal(), portSide2.ordinal());
    }
    
    /**
     * Comparator for two port contexts.
     */
    public static int comparePortContexts(final PortContext portContext1, final PortContext portContext2) {
        int portSideComparison = comparePortSides(portContext1.port.getSide(), portContext2.port.getSide());
        if (portSideComparison != 0) {
            return portSideComparison;
        }
        
        // Two ports on the same side. Depending on the side, we need an ascending or a descending comparison, because
        // ports are numbered in a clockwise manner, but we want them to be sorted left-to-right / top-to-bottom
        switch (portContext1.port.getSide()) {
        case NORTH:
        case EAST:
            return Integer.compare(portContext1.port.getVolatileId(), portContext2.port.getVolatileId());
        
        case SOUTH:
        case WEST:
            return Integer.compare(portContext2.port.getVolatileId(), portContext1.port.getVolatileId());
            
        default:
            // This should never happen since we have already checked whether a port has an undefined side in the
            // constructor of this class
            assert false;
        }
        
        return 0;
    }
    
}