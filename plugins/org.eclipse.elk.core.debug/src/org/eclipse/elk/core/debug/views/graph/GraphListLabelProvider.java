/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
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
package org.eclipse.elk.core.debug.views.graph;

import org.eclipse.elk.core.debug.ElkDebugPlugin;
import org.eclipse.elk.core.util.LoggedGraph;
import org.eclipse.jface.viewers.LabelProvider;
import org.eclipse.swt.graphics.Image;

/**
 * Label provider for graph lists.
 */
public class GraphListLabelProvider extends LabelProvider {

    /** Path for image used for elk graphs. */
    private static final String ELKGRAPH_IMAGE_PATH = "/icons/graph.gif";
    /** Path for image used for json graphs. */
    private static final String JSON_IMAGE_PATH = "/icons/graph_json.png";
    /** Path for image used for dot graphs. */
    private static final String DOT_IMAGE_PATH = "/icons/graph_dot.png";

    /** The image used for for ELK graphs. */
    private Image elkImage;
    /** The image used for for JSON graphs. */
    private Image jsonImage;
    /** The image used for for DOT graphs. */
    private Image dotImage;

    /**
     * Creates a graph list label provider.
     */
    public GraphListLabelProvider() {
        // loading icons for different log containers
        elkImage = ElkDebugPlugin.imageDescriptorFromPlugin(
                ElkDebugPlugin.PLUGIN_ID, ELKGRAPH_IMAGE_PATH).createImage();
        jsonImage = ElkDebugPlugin.imageDescriptorFromPlugin(
                ElkDebugPlugin.PLUGIN_ID, JSON_IMAGE_PATH).createImage();
        dotImage = ElkDebugPlugin.imageDescriptorFromPlugin(
                ElkDebugPlugin.PLUGIN_ID, DOT_IMAGE_PATH).createImage();
    }

    @Override
    public void dispose() {
        super.dispose();
        
        if (elkImage != null) {
            elkImage.dispose();
            elkImage = null;
        }
        
        if (jsonImage != null) {
            jsonImage.dispose();
            jsonImage = null;
        }
        
        if (dotImage != null) {
            dotImage.dispose();
            dotImage = null;
        }
    }

    @Override
    public String getText(Object element) {
        if (element instanceof LoggedGraph) {
            return ((LoggedGraph) element).getTag();            
        } else {
            return null;
        }
    }

    @Override
    public Image getImage(Object element) {
        if (element instanceof LoggedGraph) {
            switch (((LoggedGraph) element).getGraphType()) {
            case ELK:
                return elkImage;
            case JSON:
                return jsonImage;
            case DOT:
                return dotImage;
            }  
        }   
        return null;
    }
}

