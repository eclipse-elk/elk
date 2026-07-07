/*******************************************************************************
 * Copyright (c) 2008, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.ui.preferences;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.eclipse.elk.core.service.DiagramLayoutEngine;
import org.eclipse.elk.core.service.ElkServicePlugin;
import org.eclipse.elk.core.ui.ElkUiPlugin;
import org.eclipse.elk.core.ui.LayoutHandler;
import org.eclipse.elk.core.ui.Messages;
import org.eclipse.elk.core.ui.views.LayoutViewPart;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.jface.layout.LayoutConstants;
import org.eclipse.jface.preference.IPreferenceStore;
import org.eclipse.jface.preference.PreferencePage;
import org.eclipse.swt.SWT;
import org.eclipse.swt.events.SelectionAdapter;
import org.eclipse.swt.events.SelectionEvent;
import org.eclipse.swt.layout.FillLayout;
import org.eclipse.swt.layout.GridData;
import org.eclipse.swt.layout.GridLayout;
import org.eclipse.swt.program.Program;
import org.eclipse.swt.widgets.Button;
import org.eclipse.swt.widgets.Composite;
import org.eclipse.swt.widgets.Control;
import org.eclipse.swt.widgets.Group;
import org.eclipse.ui.IWorkbench;
import org.eclipse.ui.IWorkbenchPreferencePage;

/**
 * Preference page for general ELK preferences.
 */
public class LayoutPreferencePage extends PreferencePage implements IWorkbenchPreferencePage {

    private Button animationCheckBox;
    private Button zoomCheckBox;
    private Button progressCheckBox;
    private Button enableLoggingCheckBox;
    private Button storeLogsCheckBox;
    private Button execTimeCheckBox;
    private Button openLogFolderButton;

    /**
     * Creates the layout preference page.
     */
    public LayoutPreferencePage() {
        super();
    }
    
    
    /** vertical spacing between out group boxes. */
    private static final int VERTICAL_LAYOUT_SPACING = 10;

    @Override
    protected Control createContents(final Composite parent) {
        Composite composite = new Composite(parent, SWT.NONE);
        GridLayout compositeLayout = new GridLayout(1, false);
        compositeLayout.verticalSpacing = VERTICAL_LAYOUT_SPACING;
        composite.setLayout(compositeLayout);
        
        Group generalGroup = createGeneralOptionsGroup(composite);
        generalGroup.setLayoutData(new GridData(SWT.FILL, SWT.TOP, true, false));
        
        Group devOptionsGroup = createDeveloperOptionsGroup(composite);
        devOptionsGroup.setLayoutData(new GridData(SWT.FILL, SWT.TOP, true, false));
        
        return composite;
    }


    /** margin width for layouts. */
    private static final int MARGIN_WIDTH = 10;
    /** margin height for layouts. */
    private static final int MARGIN_HEIGHT = 5;
    
    /**
     * Creates the group for general options.
     * 
     * @param parent the parent control
     * @return a group with general options
     */
    private Group createGeneralOptionsGroup(final Composite parent) {
        IPreferenceStore mainPrefStore = getPreferenceStore();
        Group generalGroup = new Group(parent, SWT.NONE);
        generalGroup.setText(Messages.getString("LayoutPreferencePage.generalGroup.text")); //$NON-NLS-1$
        
        // add checkbox for animation
        animationCheckBox = new Button(generalGroup, SWT.CHECK | SWT.LEFT);
        animationCheckBox.setText(Messages.getString("LayoutPreferencePage.animations.text")); //$NON-NLS-1$
        animationCheckBox.setToolTipText(Messages.getString("LayoutPreferencePage.animations.tip")); //$NON-NLS-1$
        animationCheckBox.setSelection(mainPrefStore.getBoolean(LayoutHandler.PREF_ANIMATION));
        
        // add checkbox for zoom-to-fit
        zoomCheckBox = new Button(generalGroup, SWT.CHECK | SWT.LEFT);
        zoomCheckBox.setText(Messages.getString("LayoutPreferencePage.zoom.text")); //$NON-NLS-1$
        zoomCheckBox.setToolTipText(Messages.getString("LayoutPreferencePage.zoom.tip")); //$NON-NLS-1$
        zoomCheckBox.setSelection(mainPrefStore.getBoolean(LayoutHandler.PREF_ZOOM));
        
        // add checkbox for progress dialog
        progressCheckBox = new Button(generalGroup, SWT.CHECK | SWT.LEFT);
        progressCheckBox.setText(Messages.getString("LayoutPreferencePage.progress.text")); //$NON-NLS-1$
        progressCheckBox.setToolTipText(Messages.getString("LayoutPreferencePage.progress.tip")); //$NON-NLS-1$
        progressCheckBox.setSelection(mainPrefStore.getBoolean(LayoutHandler.PREF_PROGRESS));
        
        FillLayout layout = new FillLayout(SWT.VERTICAL);
        layout.marginWidth = MARGIN_WIDTH;
        layout.marginHeight = MARGIN_HEIGHT;
        layout.spacing = LayoutConstants.getSpacing().y;
        generalGroup.setLayout(layout);
        return generalGroup;
    }

    /**
     * Creates the group for developer options.
     * 
     * @param parent the parent control
     * @return a group with developer options
     */
    private Group createDeveloperOptionsGroup(final Composite parent) {
        IPreferenceStore servicePrefStore = ElkServicePlugin.getInstance().getPreferenceStore();
        
        Group developerGroup = new Group(parent, SWT.NONE);
        developerGroup.setText(Messages.getString("LayoutPreferencePage.developerGroup.text")); //$NON-NLS-1$
        
        GridLayout layout = new GridLayout(1, false);
        layout.marginWidth = MARGIN_WIDTH;
        layout.marginHeight = MARGIN_HEIGHT;
        layout.verticalSpacing = LayoutConstants.getSpacing().y;
        developerGroup.setLayout(layout);
        
        GridData layoutData;
        
        // add checkbox for execution time measurements
        execTimeCheckBox = new Button(developerGroup, SWT.CHECK | SWT.LEFT);
        execTimeCheckBox.setText(Messages.getString("LayoutPreferencePage.execTime.text")); //$NON-NLS-1$
        execTimeCheckBox.setToolTipText(Messages.getString("LayoutPreferencePage.execTime.tip")); //$NON-NLS-1$
        execTimeCheckBox.setSelection(servicePrefStore.getBoolean(
                DiagramLayoutEngine.PREF_DEBUG_EXEC_TIME));
        
        // add checkbox for debug graph output
        enableLoggingCheckBox = new Button(developerGroup, SWT.CHECK | SWT.LEFT);
        enableLoggingCheckBox.setText(
                Messages.getString("LayoutPreferencePage.enableLogging.text")); //$NON-NLS-1$
        enableLoggingCheckBox.setToolTipText(
                Messages.getString("LayoutPreferencePage.enableLogging.tip")); //$NON-NLS-1$
        enableLoggingCheckBox.setSelection(servicePrefStore.getBoolean(
                DiagramLayoutEngine.PREF_DEBUG_LOGGING));
        
        // add checkbox for debug graph output
        storeLogsCheckBox = new Button(developerGroup, SWT.CHECK | SWT.LEFT);
        storeLogsCheckBox.setText(Messages.getString("LayoutPreferencePage.storeLogs.text")); //$NON-NLS-1$
        storeLogsCheckBox.setToolTipText(Messages.getString("LayoutPreferencePage.storeLogs.tip")); //$NON-NLS-1$
        storeLogsCheckBox.setSelection(servicePrefStore.getBoolean(
                DiagramLayoutEngine.PREF_DEBUG_STORE));
        
        layoutData = new GridData();
        layoutData.horizontalIndent = LayoutConstants.getIndent();
        storeLogsCheckBox.setLayoutData(layoutData);
        
        // add button to open the ELK debug folder
        openLogFolderButton = new Button(developerGroup, SWT.PUSH);
        openLogFolderButton.setText(Messages.getString("LayoutPreferencePage.openLogFolder.text")); //$NON-NLS-1$
        openLogFolderButton.setToolTipText(Messages.getString("LayoutPreferencePage.openLogFolder.tip")); //$NON-NLS-1$
        openLogFolderButton.addSelectionListener(new SelectionAdapter() {
            @Override
            public void widgetSelected(final SelectionEvent event) {
                Path debugPath = Paths.get(ElkUtil.debugFolderPath());
                if (Files.isDirectory(debugPath)) {
                    Program.launch(debugPath.toString());
                }
            }
        });

        layoutData = new GridData();
        layoutData.verticalIndent = LayoutConstants.getIndent();
        openLogFolderButton.setLayoutData(layoutData);
        
        return developerGroup;
    }
    
    @Override
    public void init(final IWorkbench workbench) {
        setPreferenceStore(ElkUiPlugin.getInstance().getPreferenceStore());
    }

    @Override
    protected void performDefaults() {
        super.performDefaults();
        IPreferenceStore mainPrefStore = getPreferenceStore();
        IPreferenceStore servicePrefStore = ElkServicePlugin.getInstance().getPreferenceStore();
        
        // set default values for the general options
        animationCheckBox.setSelection(mainPrefStore.getDefaultBoolean(LayoutHandler.PREF_ANIMATION));
        zoomCheckBox.setSelection(mainPrefStore.getDefaultBoolean(LayoutHandler.PREF_ZOOM));
        progressCheckBox.setSelection(mainPrefStore.getDefaultBoolean(LayoutHandler.PREF_PROGRESS));
        execTimeCheckBox.setSelection(servicePrefStore.getDefaultBoolean(
                DiagramLayoutEngine.PREF_DEBUG_EXEC_TIME));
        enableLoggingCheckBox.setSelection(servicePrefStore.getDefaultBoolean(
                DiagramLayoutEngine.PREF_DEBUG_LOGGING));
        storeLogsCheckBox.setSelection(servicePrefStore.getDefaultBoolean(
                DiagramLayoutEngine.PREF_DEBUG_STORE));
    }

    @Override
    public boolean performOk() {
        IPreferenceStore mainPrefStore = getPreferenceStore();
        IPreferenceStore servicePrefStore = ElkServicePlugin.getInstance().getPreferenceStore();
        
        // set new values for the general options
        mainPrefStore.setValue(LayoutHandler.PREF_ANIMATION, animationCheckBox.getSelection());
        mainPrefStore.setValue(LayoutHandler.PREF_ZOOM, zoomCheckBox.getSelection());
        mainPrefStore.setValue(LayoutHandler.PREF_PROGRESS, progressCheckBox.getSelection());
        servicePrefStore.setValue(DiagramLayoutEngine.PREF_DEBUG_EXEC_TIME, execTimeCheckBox.getSelection());
        servicePrefStore.setValue(DiagramLayoutEngine.PREF_DEBUG_LOGGING, enableLoggingCheckBox.getSelection());
        servicePrefStore.setValue(DiagramLayoutEngine.PREF_DEBUG_STORE, storeLogsCheckBox.getSelection());
        
        LayoutViewPart layoutView = LayoutViewPart.findView();
        if (layoutView != null) {
            layoutView.refresh();
        }
        return true;
    }

}
