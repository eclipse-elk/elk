/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
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
package org.eclipse.elk.core.ui;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.eclipse.elk.core.data.ILayoutMetaData;
import org.eclipse.elk.core.data.LayoutAlgorithmData;
import org.eclipse.elk.core.data.LayoutCategoryData;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.util.Maybe;
import org.eclipse.jface.dialogs.Dialog;
import org.eclipse.jface.dialogs.IDialogConstants;
import org.eclipse.jface.resource.FontDescriptor;
import org.eclipse.jface.resource.ImageDescriptor;
import org.eclipse.jface.viewers.DoubleClickEvent;
import org.eclipse.jface.viewers.IDoubleClickListener;
import org.eclipse.jface.viewers.ISelectionChangedListener;
import org.eclipse.jface.viewers.ISelectionProvider;
import org.eclipse.jface.viewers.IStructuredSelection;
import org.eclipse.jface.viewers.LabelProvider;
import org.eclipse.jface.viewers.SelectionChangedEvent;
import org.eclipse.jface.viewers.StructuredSelection;
import org.eclipse.jface.viewers.TreeViewer;
import org.eclipse.jface.viewers.Viewer;
import org.eclipse.jface.viewers.ViewerComparator;
import org.eclipse.jface.viewers.ViewerFilter;
import org.eclipse.swt.SWT;
import org.eclipse.swt.events.FocusEvent;
import org.eclipse.swt.events.FocusListener;
import org.eclipse.swt.events.ModifyEvent;
import org.eclipse.swt.events.ModifyListener;
import org.eclipse.swt.graphics.Color;
import org.eclipse.swt.graphics.GC;
import org.eclipse.swt.graphics.Image;
import org.eclipse.swt.graphics.RGB;
import org.eclipse.swt.graphics.Rectangle;
import org.eclipse.swt.layout.GridData;
import org.eclipse.swt.layout.GridLayout;
import org.eclipse.swt.widgets.Button;
import org.eclipse.swt.widgets.Composite;
import org.eclipse.swt.widgets.Control;
import org.eclipse.swt.widgets.Display;
import org.eclipse.swt.widgets.Label;
import org.eclipse.swt.widgets.Shell;
import org.eclipse.swt.widgets.Text;
import org.eclipse.ui.plugin.AbstractUIPlugin;

import com.google.common.base.Strings;

/**
 * A dialog to browse and select layout algorithms or layout types.
 *
 * @author msp
 */
public class AlgorithmSelectionDialog extends Dialog {

    /** the current layouter hint as selected by the user. */
    private String layouterHint;
    /** the label for displaying the name of the current hint. */
    private Label displayNameLabel;
    /** the label for displaying the description of the current hint. */
    private Label descriptionLabel;
    /** the label for displaying the preview image. */
    private Label imageLabel;
    /** the selection provider for layout algorithms and types. */
    private ISelectionProvider selectionProvider;
    /** the cached preview images. */
    private final Map<ILayoutMetaData, Image> imageCache = new HashMap<ILayoutMetaData, Image>();
    /** the selection listeners that are added to the tree viewer when it is created. */
    private List<ISelectionChangedListener> selectionListeners
            = new LinkedList<ISelectionChangedListener>();
    
    /**
     * Creates a layout algorithm dialog.
     * 
     * @param parentShell the parent shell
     * @param currentHint the currently active layouter hint
     */
    public AlgorithmSelectionDialog(final Shell parentShell, final String currentHint) {
        super(parentShell);
        this.layouterHint = currentHint;
    }
    
    @Override
    protected void configureShell(final Shell shell) {
        super.configureShell(shell);
        shell.setText(Messages.getString("elk.ui.58")); //$NON-NLS-1$
    }
    
    /**
     * Add a selection listener that is notified when an algorithm or type is selected.
     * 
     * @param listener a listener
     */
    public void addAlgorithmSelectionListener(final ISelectionChangedListener listener) {
        if (selectionProvider != null) {
            selectionProvider.addSelectionChangedListener(listener);
        } else {
            selectionListeners.add(listener);
        }
    }
    
    /**
     * Remove a selection listener for algorithm or type changes.
     * 
     * @param listener a listener
     */
    public void removeAlgorithmSelectionListener(final ISelectionChangedListener listener) {
        if (selectionProvider != null) {
            selectionProvider.removeSelectionChangedListener(listener);
        } else {
            selectionListeners.remove(listener);
        }
    }
    
    @Override
    public boolean close() {
        imageLabel.setImage(null);
        for (Image image : imageCache.values()) {
            image.dispose();
        }
        imageCache.clear();
        return super.close();
    }
    
    /**
     * Update the currently displayed value of the description area according
     * to the tree selection.
     * 
     * @param layoutData the currently selected layout data
     */
    private void updateValue(final ILayoutMetaData layoutData) {
        layouterHint = layoutData.getId();

        String name = layoutData.getName();
        if (name == null || name.length() == 0) {
            name = layoutData instanceof LayoutAlgorithmData
                    ? Messages.getString("elk.ui.61")
                    : Messages.getString("elk.ui.8");
        }
        displayNameLabel.setText(name);
        String description = layoutData.getDescription();
        if (description == null || description.length() == 0) {
            description = Messages.getString("elk.ui.60");
        }
        descriptionLabel.setText(description);
        Image image = imageCache.get(layoutData);
        if (image == null && layoutData instanceof LayoutAlgorithmData) {
            LayoutAlgorithmData algorithmData = (LayoutAlgorithmData) layoutData;
            String path = algorithmData.getPreviewImagePath();
            String bundleId = algorithmData.getDefiningBundleId();
            if (!Strings.isNullOrEmpty(path) && !Strings.isNullOrEmpty(bundleId)) {
                ImageDescriptor imageDescriptor = AbstractUIPlugin.imageDescriptorFromPlugin(bundleId, path);
                if (imageDescriptor != null) {
                    image = createAndScaleImage(imageDescriptor);
                    if (image != null) {
                        imageCache.put(layoutData, image);
                    }
                }
            }
        }
        imageLabel.setImage(image);
        imageLabel.getParent().layout();
        
        // Enable the OK button only if a layout algorithm is selected
        Button okButton = getButton(IDialogConstants.OK_ID);
        if (okButton != null) {
            okButton.setEnabled(layoutData instanceof LayoutAlgorithmData);
        }
    }
    
    private static final RGB PREVIEW_IMG_FRAME_COLOR = new RGB(180, 180, 180);
    
    /**
     * Create an image from the image descriptor. Scale it such that it properly fits the dialog and 
     * use a {@link GC} for this; because we like anti aliasing.
     */
    private Image createAndScaleImage(final ImageDescriptor imageDescriptor) {
        // initially create an image
        Image image = imageDescriptor.createImage(false);
        if (image == null) {
            return null;
        }
        
        // scale it reasonably, and shiny
        double scale = 1;
        Rectangle bounds = image.getBounds();
        if (bounds.width > DESCRIPTION_WIDTH || bounds.height > IMAGE_MAX_HEIGHT) {
            double widthScale = Math.min(1.0, IMAGE_MAX_HEIGHT / (double) bounds.height);
            double heightScale = Math.min(1.0, DESCRIPTION_WIDTH / (double) bounds.width);
            scale = Math.min(widthScale, heightScale);
        }
        int newWidth = (int) (scale * bounds.width);
        int newHeight = (int) (scale * bounds.height);
        int imgXPadding = (DESCRIPTION_WIDTH - newWidth) / 2;
        
        Image scaled = new Image(Display.getDefault(), DESCRIPTION_WIDTH, newHeight);
        GC gc = new GC(scaled);
        gc.setAntialias(SWT.ON);
        gc.setInterpolation(SWT.HIGH);
        gc.drawImage(image, 
                0, 0, bounds.width, bounds.height, 
                imgXPadding, 0, DESCRIPTION_WIDTH - 2 * imgXPadding, newHeight - 1);

        // draw a border
        Color frameColor = new Color(Display.getCurrent(), PREVIEW_IMG_FRAME_COLOR);
        gc.setForeground(frameColor);
        gc.setAntialias(SWT.OFF);
        for (int x : new int[] {0, DESCRIPTION_WIDTH - 1}) {
            for (int y = 0; y < newHeight; ++y) {
                gc.drawPoint(x, y);
            }
        }
        for (int y : new int[] {0, newHeight - 1}) {
            for (int x = 0; x < DESCRIPTION_WIDTH; ++x) {
                gc.drawPoint(x, y);
            }
        }
        
        // clean up our mess
        gc.dispose();
        image.dispose();

        return scaled;
    }
    
    @Override
    protected Control createDialogArea(final Composite parent) {
        Composite composite = (Composite) super.createDialogArea(parent);
        ((GridLayout) composite.getLayout()).numColumns = 2;
        createSelectionTree(composite);
        createDescriptionArea(composite);    
        
        // set the initial selection according to the current layouter hint
        if (layouterHint != null) {
            LayoutMetaDataService layoutServices = LayoutMetaDataService.getInstance();
            ILayoutMetaData layoutData = layoutServices.getAlgorithmData(layouterHint);
            if (layoutData == null) {
                layoutData = layoutServices.getCategoryData(layouterHint);
            }
            if (layoutData != null) {
                selectionProvider.setSelection(new StructuredSelection(layoutData));
            }
        }
        return composite;
    }
    
    /** minimal width of the selection area. */
    private static final int SELECTION_WIDTH = 220;
    
    /**
     * Create the dialog area that displays the selection tree and filter text.
     * 
     * @param parent the parent composite
     * @return the control for the selection area
     */
    private Control createSelectionTree(final Composite parent) {
        Composite composite = new Composite(parent, SWT.NONE);
        
        // create filter text
        final Text filterText = new Text(composite, SWT.BORDER);
        filterText.setText(Messages.getString("elk.ui.59"));
        filterText.setLayoutData(new GridData(SWT.FILL, SWT.TOP, true, false));
        filterText.setForeground(parent.getDisplay().getSystemColor(SWT.COLOR_GRAY));
        
        // create tree viewer
        final TreeViewer treeViewer = new TreeViewer(composite, SWT.SINGLE | SWT.V_SCROLL | SWT.BORDER);
        final AlgorithmContentProvider contentProvider = new AlgorithmContentProvider();
        treeViewer.setContentProvider(contentProvider);
        treeViewer.setLabelProvider(new LabelProvider() {
            @Override
            public String getText(final Object element) {
                if (element instanceof LayoutAlgorithmData) {
                    LayoutAlgorithmData algoData = (LayoutAlgorithmData) element;
                    String bundleName = algoData.getBundleName();
                    if (bundleName == null) {
                        return algoData.getName();
                    } else {
                        return algoData.getName() + " (" + bundleName + ")";
                    }
                } else if (element instanceof LayoutCategoryData) {
                    LayoutCategoryData typeData = (LayoutCategoryData) element;
                    if (typeData.getName() == null) {
                        return "Other";
                    } else {
                        return typeData.getName();
                    }
                }
                return super.getText(element);
            }
        });
        treeViewer.setComparator(new ViewerComparator() {
            public int category(final Object element) {
                if (element instanceof LayoutCategoryData) {
                    LayoutCategoryData typeData = (LayoutCategoryData) element;
                    // the "Other" layout type has empty identifier and is put to the bottom
                    return typeData.getId().length() == 0 ? 1 : 0;
                }
                return super.category(element);
            }
        });
        treeViewer.getControl().setLayoutData(new GridData(SWT.FILL, SWT.FILL, true, true));
        treeViewer.setInput(LayoutMetaDataService.getInstance());
        treeViewer.expandAll();
        treeViewer.addDoubleClickListener(new IDoubleClickListener() {
            public void doubleClick(final DoubleClickEvent event) {
                okPressed();
            }
        });
        
        // set up a filter on the tree viewer using the filter text
        final Maybe<Boolean> filterChanged = new Maybe<Boolean>(Boolean.FALSE);
        final Maybe<Boolean> filterLeft = new Maybe<Boolean>(Boolean.FALSE);
        filterText.addModifyListener(new ModifyListener() {
            public void modifyText(final ModifyEvent e) {
                if (!filterChanged.get()) {
                    filterChanged.set(Boolean.TRUE);
                    filterText.setForeground(null);
                    int pos = filterText.getCaretPosition();
                    String newText = filterText.getText(pos - 1, pos - 1);
                    filterText.setText(newText);
                    filterText.setSelection(pos);
                } else {
                    contentProvider.updateFilter(filterText.getText());
                    treeViewer.refresh();
                    treeViewer.expandAll();
                    ILayoutMetaData selected = contentProvider.getBestFilterMatch();
                    if (selected != null) {
                        treeViewer.setSelection(new StructuredSelection(selected));
                    }
                }
            }
        });
        filterText.addFocusListener(new FocusListener() {
            public void focusGained(final FocusEvent e) {
                if (filterLeft.get() && !filterChanged.get()) {
                    filterChanged.set(Boolean.TRUE);
                    filterText.setForeground(null);
                    filterText.setText("");
                }
            }
            public void focusLost(final FocusEvent e) {
                filterLeft.set(Boolean.TRUE);
            }
        });
        treeViewer.addFilter(new ViewerFilter() {
            public boolean select(final Viewer viewer, final Object parentElement,
                    final Object element) {
                return contentProvider.applyFilter(element);
            }
        });
        
        // add a selection listener to the tree so that the selected element is displayed
        treeViewer.addSelectionChangedListener(new ISelectionChangedListener() {
            public void selectionChanged(final SelectionChangedEvent event) {
                IStructuredSelection selection = (IStructuredSelection) event.getSelection();
                Object element = selection.getFirstElement();
                if (element instanceof ILayoutMetaData) {
                    updateValue((ILayoutMetaData) element);
                }
            }
        });
        
        composite.setLayout(new GridLayout());
        GridData gridData = new GridData(SWT.FILL, SWT.FILL, true, true);
        gridData.minimumWidth = SELECTION_WIDTH;
        composite.setLayoutData(gridData);
        
        // register all selection listeners that have been stored in a list
        selectionProvider = treeViewer;
        for (ISelectionChangedListener listener : selectionListeners) {
            selectionProvider.addSelectionChangedListener(listener);
        }
        return composite;
    }
    
    /** width of the description area. */
    private static final int DESCRIPTION_WIDTH = 300;
    /** vertical spacing in the description area. */
    private static final int DESCR_SPACING = 12;
    /** maximum height of preview images .*/
    private static final int IMAGE_MAX_HEIGHT = 200;
    
    /**
     * Create the dialog area that displays the description of a layout algorithm.
     * 
     * @param parent the parent composite
     * @return the control for the description area
     */
    private Control createDescriptionArea(final Composite parent) {
        Composite composite = new Composite(parent, SWT.NONE);
        
        // create label for the display name
        displayNameLabel = new Label(composite, SWT.NONE);
        FontDescriptor fontDescriptor = FontDescriptor.createFrom(parent.getFont());
        fontDescriptor = fontDescriptor.increaseHeight(2).setStyle(SWT.BOLD);
        displayNameLabel.setFont(fontDescriptor.createFont(parent.getDisplay()));
        displayNameLabel.setLayoutData(new GridData(SWT.FILL, SWT.TOP, true, false));
        
        // create label for the description
        descriptionLabel = new Label(composite, SWT.WRAP);
        GridData descriptionLayoutData = new GridData(SWT.FILL, SWT.FILL, true, true);
        descriptionLayoutData.widthHint = DESCRIPTION_WIDTH;
        descriptionLabel.setLayoutData(descriptionLayoutData);
        
        // create label for the preview image
        imageLabel = new Label(composite, SWT.NONE);
        GridData imageLayoutData = new GridData(SWT.CENTER, SWT.BOTTOM, true, false);
        imageLabel.setLayoutData(imageLayoutData);
        
        GridLayout compositeLayout = new GridLayout();
        compositeLayout.verticalSpacing = DESCR_SPACING;
        composite.setLayout(compositeLayout);
        composite.setLayoutData(new GridData(SWT.FILL, SWT.FILL, true, true));
        return composite;
    }
    
    /**
     * The layouter hint that was selected by the user.
     * 
     * @return the selected layouter hint
     */
    public String getSelectedHint() {
        return layouterHint;
    }

}
