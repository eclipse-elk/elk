/*******************************************************************************
 * Copyright (c) 2019, 2020 Kiel University and others.
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
package org.eclipse.elk.core.debug.wizard;

import java.io.ByteArrayInputStream;
import java.util.Arrays;

import org.eclipse.core.resources.IFile;
import org.eclipse.core.resources.IFolder;
import org.eclipse.core.resources.IProject;
import org.eclipse.core.resources.IProjectDescription;
import org.eclipse.core.runtime.CoreException;
import org.eclipse.core.runtime.IProgressMonitor;
import org.eclipse.elk.core.debug.wizard.templates.LayoutProviderTemplate;
import org.eclipse.elk.core.debug.wizard.templates.MelkTemplate;
import org.eclipse.jface.wizard.Wizard;
import org.eclipse.pde.core.plugin.IPluginModelBase;
import org.eclipse.pde.core.plugin.IPluginReference;
import org.eclipse.pde.ui.IFieldData;
import org.eclipse.pde.ui.IPluginContentWizard;
import org.eclipse.pde.ui.templates.PluginReference;


/**
 * Plug-in project wizard for ELK Layout Algorithms.
 */
public class AlgorithmPluginContentWizard extends Wizard implements IPluginContentWizard {
    
    /** IDs of the projects all layout algorithm projects need to depend on. */
    private static final String[] PLUGIN_DEPENDENCIES = {
        "org.eclipse.elk.core",
        "org.eclipse.elk.graph",
        "org.eclipse.elk.alg.common",
        "com.google.guava"
    };
    /** ID of the Xtext nature. */
    private static final String XTEXT_NATURE = "org.eclipse.xtext.ui.shared.xtextNature";

    /** ID of the plug-in project to be created by the wizard. */
    private String id;
    /** Name of the plug-in project to be created by the wizard. */
    private String plugInName;
    /** Name of our new layout algorithm to be created in the new plug-in. */
    private String algorithmName;
    
    /* GUI */
    private AlgorithmProjectCreationPage page;


    @Override
    public void init(IFieldData data) {
        id = data.getId();
        plugInName = data.getName();
    }


    @Override
    public void addPages() {
        super.addPages();
        
        // Simply add our single page here
        page = new AlgorithmProjectCreationPage(plugInName);
        addPage(page);
    }

    @Override
    public IPluginReference[] getDependencies(String schemaVersion) {
        return Arrays.stream(PLUGIN_DEPENDENCIES)
            .map(id -> new PluginReference(id))
            .toArray(IPluginReference[]::new);
    }

    @Override
    public String[] getNewFiles() {
        // We don't want anything to be added to build.properties, so return an empty string array (returning null
        // causes an NPE)
        return new String[0];
    }

    @Override
    public boolean performFinish() {
        // We always allow the user to finish
        return page != null && page.getAlgorithmName().trim().length() > 0;
    }

    @Override
    public boolean performFinish(IProject project, IPluginModelBase model, IProgressMonitor monitor) {
        // Executes file generation. Creates packages, the melk file, the layout provider and the plugin.xml.
        
        // Get algorithm name and capitalize it
        if (page != null && !page.getAlgorithmName().trim().equals("")) {
            algorithmName = page.getAlgorithmName().substring(0, 1).toUpperCase() 
                    + page.getAlgorithmName().substring(1);
        } else {
            algorithmName = plugInName.substring(0, 1).toUpperCase() + plugInName.substring(1);
        }

        // Create the base package
        String path = "src/";
        try {
            for (String part : id.split("\\.")) {
                IFolder folder = project.getFolder(path + part);
                folder.create(false, false, monitor);
                path = path + part + "/";
            }
        } catch (CoreException e) {
            e.printStackTrace();
            return false;
        }

        // Populate the project
        createLayoutProvider(project, monitor);
        
        createMelkFile(project, monitor);
        createServiceFile(project, monitor);
        addXtextNature(project, monitor);
        
        return true;
    }

    /**
     * Creates the layout provider Java class from template.
     */
    private void createLayoutProvider(IProject project, IProgressMonitor monitor) {
        String fileContent = LayoutProviderTemplate.buildFileContent(
                id,
                algorithmName);

        String path = "/src/" + id.replace(".", "/") + "/" + algorithmName + "LayoutProvider.java";
        IFile layoutProvider = project.getFile(path);
        try {
            layoutProvider.create(new ByteArrayInputStream(fileContent.getBytes()), false, monitor);
        } catch (CoreException e) {
            e.printStackTrace();
        }
    }

    /**
     * Creates the project's MELK file from template.
     */
    private void createMelkFile(IProject project, IProgressMonitor monitor) {
        String fileContent = MelkTemplate.buildFileContent(
                id,
                id,
                algorithmName + "LayoutProvider",
                id,
                algorithmName,
                algorithmName);

        IFile melk = project.getFile("/src/" + id.replace(".", "/") + "/" + plugInName + ".melk");
        try {
            melk.create(new ByteArrayInputStream(fileContent.getBytes()), false, monitor);
        } catch (CoreException e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Add layout provider extension to plug-in to register algorithm with ELK.
     */
    private void createServiceFile(IProject project, IProgressMonitor monitor) {
        String fileContent = id + "." + algorithmName + "MetadataProvider";
        
        // Create necessary folders
        try {
            IFolder srcMetaInf = project.getFolder("/src/META-INF");
            srcMetaInf.create(true, true, monitor);
            
            IFolder srcMetaInfServices = project.getFolder("/src/META-INF/services");
            srcMetaInfServices.create(true, true, monitor);
            
            IFile serviceFile = project.getFile(
                    "/src/META-INF/services/org.eclipse.elk.core.data.ILayoutMetaDataProvider");
            serviceFile.create(new ByteArrayInputStream(fileContent.getBytes()), false, monitor);
        } catch (CoreException e) {
            e.printStackTrace();
        }
    }

    /**
     * Adds the Xtext nature to the project to ensure that MELK files are compiled.
     */
    private void addXtextNature(IProject project, IProgressMonitor monitor) {
        try {
            IProjectDescription description = project.getDescription();
            String[] natures = description.getNatureIds();

            String[] newNatures = new String[natures.length + 1];
            System.arraycopy(natures, 0, newNatures, 0, natures.length);
            newNatures[natures.length] = XTEXT_NATURE;

            description.setNatureIds(newNatures);
            project.setDescription(description, monitor);
        } catch (CoreException e) {
            e.printStackTrace();
        }
    }
}
