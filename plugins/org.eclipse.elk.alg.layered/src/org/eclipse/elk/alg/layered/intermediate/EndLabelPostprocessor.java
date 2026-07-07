/*******************************************************************************
 * Copyright (c) 2012, 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import java.util.Map;

import org.eclipse.elk.alg.common.nodespacing.cellsystem.LabelCell;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * <p>After the {@link EndLabelPreprocessor} has done all the major work for us, each node may have a list of label
 * cells full of edge end labels associated with it, along with proper label cell coordinates. The problem is that
 * when the preprocessor was run, the final coordinates of all the nodes (and with that, the final coordinates of all
 * the edges) were not known yet. Now they are, so the label cell coordinates can be offset accordingly and label
 * coordinates can be applied.</p>
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>a layered graph</dd>
 *     <dd>end labels of all edges are contained in label cells associated with nodes</dd>
 *     <dd>label cells have correct sizes and settings and are placed relative to their node's upper left corner</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>end labels have proper coordinates assigned to them</dd>
 *   <dt>Slots:</dt>
 *     <dd>After phase 5.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>{@link EndLabelSorter}</dd>
 * </dl>
 */
public final class EndLabelPostprocessor implements ILayoutProcessor<LGraph> {
    
    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("End label post-processing", 1);
        
        // We iterate over each node's label cells and offset and place them
        layeredGraph.getLayers().stream()
                .flatMap(layer -> layer.getNodes().stream())
                .filter(node -> ((node.getType() == NodeType.NORMAL || node.getType() == NodeType.EXTERNAL_PORT) && node.hasProperty(InternalProperties.END_LABELS)))
                .forEach(node -> processNode(node));
        
        monitor.done();
    }
    
    private void processNode(final LNode node) {
        assert node.hasProperty(InternalProperties.END_LABELS);
        
        // The node should have a non-empty list of label cells, or something went TERRIBLY WRONG!!!
        Map<LPort, LabelCell> endLabelCells = node.getProperty(InternalProperties.END_LABELS);
        assert !endLabelCells.isEmpty();
        
        KVector nodePos = node.getPosition();
        
        for (LabelCell labelCell : endLabelCells.values()) {
            ElkRectangle labelCellRect = labelCell.getCellRectangle();
            labelCellRect.move(nodePos);
            
            labelCell.applyLabelLayout();
        }
        
        // Remove label cells
        node.setProperty(InternalProperties.END_LABELS, null);
    }

}
