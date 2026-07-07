/*******************************************************************************
 * Copyright (c) 2018 TypeFox GmbH and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.text.ui.layout;

import org.eclipse.elk.core.service.IDiagramLayoutConnector;
import org.eclipse.elk.core.service.LayoutMapping;
import org.eclipse.elk.core.ui.rendering.GraphRenderingConfigurator;
import org.eclipse.elk.core.ui.rendering.LayoutGraphDialog;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.properties.IPropertyHolder;
import org.eclipse.emf.ecore.EObject;
import org.eclipse.emf.ecore.util.EcoreUtil;
import org.eclipse.swt.widgets.Display;
import org.eclipse.ui.IWorkbenchPart;
import org.eclipse.ui.PlatformUI;
import org.eclipse.xtext.ui.editor.XtextEditor;

/**
 * This connector copies the graph from an Xtext editor and displays the layout result in a dialog.
 */
public class TextLayoutConnector implements IDiagramLayoutConnector {

    @Override
    public LayoutMapping buildLayoutGraph(IWorkbenchPart workbenchPart, Object diagramPart) {
        if (workbenchPart instanceof XtextEditor) {
            XtextEditor xtextEditor = (XtextEditor) workbenchPart;
            return xtextEditor.getDocument().readOnly(resource -> {
                EObject content = resource.getContents().isEmpty() ? null : resource.getContents().get(0);
                if (content instanceof ElkNode) {
                    return buildLayoutGraph(xtextEditor, (ElkNode) content);
                } else {
                    throw new IllegalArgumentException("Not supported by this layout connector: " + content);
                }
            });
        } else {
            throw new IllegalArgumentException("Not supported by this layout connector: " + workbenchPart);
        }
    }
    
    protected LayoutMapping buildLayoutGraph(XtextEditor xtextEditor, ElkNode rootNode) {
        LayoutMapping mapping = new LayoutMapping(xtextEditor);
        mapping.setParentElement(rootNode);
        mapping.setLayoutGraph(EcoreUtil.copy(rootNode));
        return mapping;
    }

    @Override
    public void applyLayout(LayoutMapping mapping, IPropertyHolder settings) {
        Display display = PlatformUI.getWorkbench().getDisplay();
        display.asyncExec(() -> {
            EObject parentElement = (EObject) mapping.getParentElement();
            GraphRenderingConfigurator configurator = new GraphRenderingConfigurator(display);
            LayoutGraphDialog dialog = new LayoutGraphDialog(mapping.getLayoutGraph(),
                    parentElement.eResource(), configurator, display.getActiveShell());
            dialog.open();
        });
    }

}
