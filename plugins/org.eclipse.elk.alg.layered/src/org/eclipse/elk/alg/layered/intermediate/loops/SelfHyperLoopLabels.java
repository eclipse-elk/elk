/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.loops;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.eclipse.elk.alg.layered.graph.LGraphUtil;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.intermediate.LabelManagementProcessor;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.labels.ILabelManager;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.options.PortSide;

/**
 * Represents the labels associated with a {@link SelfHyperLoop}. The labels are gathered as edges are added to the
 * hyper loop. Their actual position is determined by a port side and an alignment on that side.
 */
public class SelfHyperLoopLabels {
    
    /**
     * Describes all the different alignments a label might have. Most are alignments for labels on the northern or
     * southern side. Labels on the eastern or western side are restricted to {@link Alignment#TOP}. 
     */
    public static enum Alignment {
        /** A northern or southern centered label. */
        CENTER,
        /** A northern or southern left-aligned label. */
        LEFT,
        /** A northern or southern right-aligned label. */
        RIGHT,
        /** An eastern or western top-aligned label. */
        TOP
    }

    // SUPPRESS CHECKSTYLE NEXT 2 VisibilityModifier
    /** An ID that can be used arbitrarily. */
    public int id = 0;
    
    /** The labels represented by this instance. */
    private final List<LLabel> lLabels = new ArrayList<>();
    /** The size required to place all labels. */
    private final KVector size = new KVector();
    /** The position of our bounding box's top left corner. */
    private final KVector position = new KVector();
    /** The graph's layout direction. Influences the way the labels will be placed. */
    private final Direction layoutDirection;
    /** Space to leave between adjacent labels. */
    private final double labelLabelSpacing;
    
    /** The side the label is placed on. */
    private PortSide side;
    /** The label's alignment. */
    private Alignment alignment;
    /** The port any non-center alignment is relative to. */
    private SelfLoopPort alignmentReferenceSLPort;
    
    
    /**
     * Creates a new instance for the given {@link SelfHyperLoop}.
     */
    public SelfHyperLoopLabels(final SelfHyperLoop slLoop) {
        // Initialize our properties from the graph
        LNode lNode = slLoop.getSLHolder().getLNode();
        
        layoutDirection = lNode.getGraph().getProperty(LayeredOptions.DIRECTION);
        labelLabelSpacing = LGraphUtil.getIndividualOrInherited(lNode, LayeredOptions.SPACING_LABEL_LABEL);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // LLabel Access
    
    /**
     * Adds the given collection of labels to this instance.
     */
    public void addLLabels(final Collection<LLabel> newLLabels) {
        for (LLabel newLLabel : newLLabels) {
            lLabels.add(newLLabel);
            updateSize(newLLabel);
        }
    }

    /**
     * Returns the labels represented by this instance.
     */
    public List<LLabel> getLLabels() {
        return lLabels;
    }
    
    /**
     * Updates our size under the assumption that the given label was just added to this instance.
     */
    private void updateSize(final LLabel newLLabel) {
        assert lLabels.get(lLabels.size() - 1) == newLLabel;
        
        KVector newLLabelSize = newLLabel.getSize();
        
        if (layoutDirection.isHorizontal()) {
            // The labels will be stacked vertically
            size.x = Math.max(size.x, newLLabelSize.x);
            size.y += newLLabelSize.y;
            
            // Add a label-label spacing if we already had labels
            if (lLabels.size() > 1) {
                size.y += labelLabelSpacing;
            }
            
        } else {
            // The labels will be stacked horizontally
            size.x += newLLabelSize.x;
            size.y = Math.max(size.y, newLLabelSize.y);
            
            // Add a label-label spacing if we already had labels
            if (lLabels.size() > 1) {
                size.x += labelLabelSpacing;
            }
        }
    }
    
    /**
     * Applies label management with the given target width to our labels and updates our size. 
     */
    public void applyLabelManagement(final ILabelManager labelManager, final double targetWidth) {
        size.set(LabelManagementProcessor.doManageLabels(
                labelManager, lLabels, targetWidth, labelLabelSpacing, layoutDirection.isVertical()));
    }
    
    /**
     * Applies the bounding box placement to the individual labels.
     */
    public void applyPlacement(final KVector offset) {
        if (layoutDirection.isHorizontal()) {
            applyPlacementForHorizontalLayout(offset);
        } else {
            applyPlacementVerticalForVerticalLayout(offset);
        }
    }

    private void applyPlacementForHorizontalLayout(final KVector offset) {
        double x = position.x;
        double y = position.y;
        
        for (LLabel lLabel : lLabels) {
            KVector labelPos = lLabel.getPosition();
            
            // X coordinate depends on alignment and / or side
            if (alignment == Alignment.LEFT || side == PortSide.EAST) {
                labelPos.x = x;
            } else if (alignment == Alignment.RIGHT || side == PortSide.WEST) {
                labelPos.x = x + size.x - lLabel.getSize().x;
            } else {
                // Alignment is center
                labelPos.x = x + (size.x - lLabel.getSize().x) / 2;
            }
            
            labelPos.y = y;
            labelPos.add(offset);
            
            y += lLabel.getSize().y + labelLabelSpacing;
        }
    }

    private void applyPlacementVerticalForVerticalLayout(final KVector offset) {
        double x = position.x;
        double y = position.y;
        
        for (LLabel lLabel : lLabels) {
            KVector labelPos = lLabel.getPosition();
            
            labelPos.x = x;
            
            // We always top-align, except for the northern side
            if (side == PortSide.NORTH) {
                labelPos.y = y + size.y - lLabel.getSize().y;
            } else {
                labelPos.y = y;
            }
            
            labelPos.add(offset);
            
            x += lLabel.getSize().x + labelLabelSpacing;
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Label Placement
    
    /**
     * The size of the box required to place the labels. Shouldn't be modified by clients.
     */
    public KVector getSize() {
        return size;
    }
    
    /**
     * Position of the bounding box's top left corner. To be modified by clients.
     */
    public KVector getPosition() {
        return position;
    }
    
    /**
     * Returns the side this label is placed on. May be {@code null} if the side was not determined yet.
     */
    public PortSide getSide() {
        return side;
    }
    
    /**
     * Sets the side this label is placed on.
     */
    public void setSide(final PortSide side) {
        this.side = side;
    }
    
    /**
     * Returns this label's alignment May be {@code null} if the alignment was not determined yet.
     */
    public Alignment getAlignment() {
        return alignment;
    }
    
    /**
     * Sets this label's alignment.
     */
    public void setAlignment(final Alignment alignment) {
        this.alignment = alignment;
    }
    
    /**
     * Returns the port the alignment is relative to. This will only be non-{@code null} (and thus meaningful) if the
     * alignment is not {@link Alignment#CENTER}.
     */
    public SelfLoopPort getAlignmentReferenceSLPort() {
        return alignmentReferenceSLPort;
    }
    
    /**
     * Sets the port the alignment is relative to. Should only be called if the alignment is not
     * {@link Alignment#CENTER}.
     */
    public void setAlignmentReferenceSLPort(final SelfLoopPort alignmentReferenceSLPort) {
        this.alignmentReferenceSLPort = alignmentReferenceSLPort;
    }
    
}
