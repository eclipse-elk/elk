/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.debug.wizard;

import org.eclipse.jface.wizard.WizardPage;
import org.eclipse.swt.SWT;
import org.eclipse.swt.events.KeyAdapter;
import org.eclipse.swt.events.KeyEvent;
import org.eclipse.swt.layout.GridData;
import org.eclipse.swt.layout.GridLayout;
import org.eclipse.swt.widgets.Composite;
import org.eclipse.swt.widgets.Label;
import org.eclipse.swt.widgets.Text;

/**
 * Wizard Page for Layout Algorithm Project wizard. One is able to enter a layout algorithm name which can be retrieved
 * by the wizard for further use.
 */
public class AlgorithmProjectCreationPage extends WizardPage {
    
    /** Name of the plug-in to be created by the wizard. This will be the base for the layout algorithm's name. */
    private String plugInName;
    
    /* GUI */
    private Text algorithmName;
    private Composite container;

    
    /**
     * Wizard page constructor.
     * 
     * @param plugInName the plug-ins name.
     */
    public AlgorithmProjectCreationPage(String plugInName) {
        super("Layout Algorithm");
        this.plugInName = plugInName;
        setTitle("Layout Algorithm Plug-In");
        setDescription("Choose the name of the algorithm in your project");
    }

    /**
     * Sets up the view.
     */
    @Override
    public void createControl(Composite parent) {
        // general layout of the page
        container = new Composite(parent, SWT.NONE);
        GridLayout layout = new GridLayout();
        container.setLayout(layout);
        layout.numColumns = 2;
        
        // add label in first column
        Label algNameLabel = new Label(container, SWT.NONE);
        algNameLabel.setText("Algorithm name:");
        
        // add a text input in second column
        algorithmName = new Text(container, SWT.BORDER | SWT.SINGLE);
        algorithmName.setText(plugInName);
        GridData gd = new GridData(GridData.FILL_HORIZONTAL);
        algorithmName.setLayoutData(gd);
        algorithmName.addKeyListener(new KeyAdapter() {
            @Override
            public void keyPressed(KeyEvent e) {
                setPageComplete(!algorithmName.getText().isEmpty());
            }
        });
        
        setControl(container);
        setPageComplete(true);
    }
    
    /**
     * Return the name of the algorithm to be created.
     */
    public String getAlgorithmName() {
        return algorithmName.getText();
    }
}
