/*******************************************************************************
 * Copyright (c) 2008, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.graphviz.layouter.preferences;

import org.eclipse.elk.alg.graphviz.layouter.GraphvizLayoutProvider;
import org.eclipse.elk.alg.graphviz.layouter.GraphvizTool;
import org.eclipse.elk.core.data.LayoutAlgorithmData;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.jface.preference.BooleanFieldEditor;
import org.eclipse.jface.preference.FieldEditorPreferencePage;
import org.eclipse.jface.preference.FileFieldEditor;
import org.eclipse.jface.preference.IntegerFieldEditor;
import org.eclipse.swt.widgets.Composite;
import org.eclipse.ui.IWorkbench;
import org.eclipse.ui.IWorkbenchPreferencePage;

/**
 * The Graphviz preference page.
 * 
 * @author ars
 * @author msp
 */
public class GraphvizPreferencePage extends FieldEditorPreferencePage implements
        IWorkbenchPreferencePage {

    /** identifier of the preference page. */
    public static final String ID = "org.eclipse.elk.alg.graphviz.preferences";

    /**
     * Creates a Graphviz preference page.
     */
    public GraphvizPreferencePage() {
        super(GRID);
        setDescription(
                  "Controls how the Eclipse Layout Kernel interacts with the Graphviz layout tools. "
                + "The Graphviz layout tools are available at http://www.graphviz.org/. "
                + "If the 'dot' executable cannot be found in default locations, "
                + "its path must be entered here.");
    }

    @Override
    public void createFieldEditors() {
        // Process group
        Composite executableEditorParent = getFieldEditorParent();
        FileFieldEditor executableEditor = new FileFieldEditor(GraphvizTool.PREF_GRAPHVIZ_EXECUTABLE,
                "Path to 'dot' executable:", executableEditorParent);
        executableEditor.setValidateStrategy(FileFieldEditor.VALIDATE_ON_KEY_STROKE);
        addField(executableEditor);

        Composite timeoutEditorParent = getFieldEditorParent();
        IntegerFieldEditor timeoutEditor = new IntegerFieldEditor(GraphvizTool.PREF_TIMEOUT,
                "Timeout (ms):", timeoutEditorParent);
        timeoutEditor.setValidRange(GraphvizTool.PROCESS_MIN_TIMEOUT, Integer.MAX_VALUE);
        addField(timeoutEditor);

        BooleanFieldEditor restartGraphvizProcessCheckbox = new BooleanFieldEditor(
                GraphvizLayoutProvider.PREF_GRAPHVIZ_REUSE_PROCESS,
                "Reuse single graphviz process for better performance", getFieldEditorParent());
        addField(restartGraphvizProcessCheckbox);
    }

    @Override
    public void init(final IWorkbench workbench) {
        setPreferenceStore(GraphvizLayouterPreferenceStore.getInstance().getPreferenceStore());
    }
    
    @Override
    public boolean performOk() {
        // dispose all cached Graphviz instances to ensure creation of new processes
        for (LayoutAlgorithmData data : LayoutMetaDataService.getInstance().getAlgorithmData()) {
            if ("Graphviz".equals(data.getBundleName())) {
                data.getInstancePool().clear();
            }
        }
        return super.performOk();
    }

}
