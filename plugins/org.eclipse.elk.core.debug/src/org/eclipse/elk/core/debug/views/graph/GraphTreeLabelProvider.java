/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.debug.views.graph;

import java.time.format.DateTimeFormatter;

import org.eclipse.elk.core.debug.ElkDebugPlugin;
import org.eclipse.elk.core.debug.model.ExecutionInfo;
import org.eclipse.jface.viewers.LabelProvider;
import org.eclipse.jface.viewers.StyledString;
import org.eclipse.jface.viewers.DelegatingStyledCellLabelProvider.IStyledLabelProvider;
import org.eclipse.swt.graphics.Image;

/**
 * Label provider for graph containers.
 */
public class GraphTreeLabelProvider extends LabelProvider implements IStyledLabelProvider {

    /** Path for image used for containers with graphs. */
    private static final String GRAPH_IMAGE_PATH = "/icons/log.png";
    /** Path for image used for containers without graphs. */
    private static final String NO_GRAPH_IMAGE_PATH = "/icons/no_log.png";

    /** The image used for containers with graphs. */
    private Image graphImage;
    /** The image used for containers without graphs. */
    private Image noGraphImage;

    public GraphTreeLabelProvider() {
        // loading icons for different graph containers
        graphImage =
                ElkDebugPlugin.imageDescriptorFromPlugin(ElkDebugPlugin.PLUGIN_ID, GRAPH_IMAGE_PATH).createImage();
        noGraphImage =
                ElkDebugPlugin.imageDescriptorFromPlugin(ElkDebugPlugin.PLUGIN_ID, NO_GRAPH_IMAGE_PATH).createImage();
    }

    @Override
    public void dispose() {
        super.dispose();

        if (noGraphImage != null) {
            noGraphImage.dispose();
            noGraphImage = null;
        }

        if (graphImage != null) {
            graphImage.dispose();
            graphImage = null;
        }
    }

    @Override
    public StyledString getStyledText(Object element) {
        if (element instanceof ExecutionInfo) {
            ExecutionInfo execution = (ExecutionInfo) element;
            StyledString text = new StyledString(execution.getName());
            
            if (execution.getParent() == null) {
                // Add the creation time for top-level elements
                text.append(" (" + DateTimeFormatter.ISO_INSTANT.format(execution.getCreationTime()) + ")",
                        StyledString.COUNTER_STYLER);
            }
            
            return text;
            
        } else {
            return null;
        }
    }

    @Override
    public Image getImage(Object element) {
        if (element instanceof ExecutionInfo) {
            return ((ExecutionInfo) element).hasLoggedGraphs() ? graphImage : noGraphImage;
        } else {
            return null;
        }
    }

}
