/*******************************************************************************
 * Copyright (c) 2018 TypeFox GmbH and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.ui.rendering;

import org.eclipse.elk.graph.ElkNode;
import org.eclipse.emf.ecore.resource.Resource;
import org.eclipse.jface.dialogs.Dialog;
import org.eclipse.swt.layout.FillLayout;
import org.eclipse.swt.widgets.Composite;
import org.eclipse.swt.widgets.Control;
import org.eclipse.swt.widgets.Shell;

/**
 * A dialog that renders an ElkGraph.
 */
public class LayoutGraphDialog extends Dialog {
    
    private final ElkNode graph;
    private final Resource resource;
    private final GraphRenderingConfigurator configurator;
    
    private GraphRenderingCanvas graphCanvas;

    /**
     * Create a layout graph dialog.
     */
    public LayoutGraphDialog(final ElkNode graph, final Resource resource,
            final GraphRenderingConfigurator configurator, final Shell parentShell) {
        super(parentShell);
        this.graph = graph;
        this.resource = resource;
        this.configurator = configurator;
    }
    
    @Override
    protected boolean isResizable() {
        return true;
    }
    
    @Override
    protected void configureShell(final Shell shell) {
        super.configureShell(shell);
        String fileName = null;
        if (resource != null && resource.getURI() != null) {
            fileName = resource.getURI().lastSegment();
        }
        if (fileName == null || fileName.isEmpty()) {
            shell.setText("Layout Graph");
        } else {
            shell.setText(fileName);
        }
    }
    
    @Override
    protected Control createContents(final Composite parent) {
        parent.setLayout(new FillLayout());
        graphCanvas = new GraphRenderingCanvas(parent, configurator);
        graphCanvas.setupMouseInteraction();
        graphCanvas.setLayoutGraph(graph);
        return graphCanvas;
    }

}
